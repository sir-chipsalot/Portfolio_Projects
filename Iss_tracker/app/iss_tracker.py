from flask import Flask, request, jsonify
import redis
import time
from datetime import datetime
import math
import json
from geopy.geocoders import Nominatim
from astropy import coordinates, units
from astropy.time import Time
from fetch_iss_tracker import fetch_iss_data
import logging

# Set up logging
logging.basicConfig(level='DEBUG')

app = Flask(__name__)

# Redis connection
def get_redis_client():
    return redis.Redis(host='redis-db', port=6379, db=0)

rd = get_redis_client()

# Load ISS data into Redis if not already loaded
def load_iss_data():
    if not rd.exists("iss_data"):
        # Fetch data from ISS API or local file
        fetch_iss_data()
        with open("Iss_data.json", "r") as e:
            iss_data = json.load(e)
            rd.set("iss_data", json.dumps(iss_data))
            print("ISS data loaded into Redis.")
    else:
        print("Failed to fetch ISS data.")

# Helper function to calculate instantaneous speed
def instantaneous_speed(state_vector):
    try:
        x_dot = float(state_vector["X_DOT"]["#text"])
        y_dot = float(state_vector["Y_DOT"]["#text"])
        z_dot = float(state_vector["Z_DOT"]["#text"])
        return math.sqrt(x_dot**2 + y_dot**2 + z_dot**2)
    except Exception as e:
        logging.error(f"Failed to calculate speed: {e}")
        return None

# Helper function to compute latitude, longitude, and altitude
def compute_location_astropy(state_vector):
    try:
        x = float(state_vector["X"]["#text"])
        y = float(state_vector["Y"]["#text"])
        z = float(state_vector["Z"]["#text"])

        # Convert epoch to a timestamp
        epoch = state_vector["EPOCH"]
        this_epoch = time.strftime('%Y-%m-%d %H:%M:%S', time.strptime(epoch, '%Y-%jT%H:%M:%S'))

        # Compute location using astropy
        cartrep = coordinates.CartesianRepresentation([x, y, z], unit=units.km)
        gcrs = coordinates.GCRS(cartrep, obstime=this_epoch)
        itrs = gcrs.transform_to(coordinates.ITRS(obstime=this_epoch))
        loc = coordinates.EarthLocation(*itrs.cartesian.xyz)

        return loc.lat.value, loc.lon.value, loc.height.value
    except Exception as e:
        logging.error(f"Failed to compute location: {e}")
        return None, None, None

# Helper function to find the closest epoch to the current time
def find_closest_epoch(state_vectors):
    current_time = time.time()
    min_difference = float("inf")
    closest_entry = None

    for item in state_vectors:
        try:
            epoch_time = datetime.strptime(item["EPOCH"], "%Y-%jT%H:%M:%S.%fZ")
            epoch_iter = time.mktime(epoch_time.timetuple())
            time_diff = abs(epoch_iter - current_time)

            if time_diff < min_difference:
                min_difference = time_diff
                closest_entry = item
        except Exception as e:
            logging.error(f"Failed to parse epoch time: {e}")
            exit(1)

    return closest_entry

# Route: /epochs
@app.route('/epochs', methods=['GET'])
def get_epochs():
    try:
        iss_data = json.loads(rd.get("iss_data"))
        state_vectors = iss_data["ndm"]["oem"]["body"]["segment"]["data"]["stateVector"]

        # Handle query parameters: limit and offset
        limit = request.args.get('limit', default=None, type=int)
        offset = request.args.get('offset', default=0, type=int)

        if limit is not None:
            state_vectors = state_vectors[offset:offset + limit]
        else:
            state_vectors = state_vectors[offset:]

        return jsonify({"state_vectors": state_vectors})
    except Exception as e:
        logging.error(f"Failed to fetch epochs: {e}")
        return jsonify({"error": "An error occurred"})

# Route: /epochs/<epoch>/speed
@app.route('/epochs/<epoch>/speed', methods=['GET'])
def get_epoch_speed(epoch):
    try:
        iss_data = json.loads(rd.get("iss_data"))
        state_vectors = iss_data["ndm"]["oem"]["body"]["segment"]["data"]["stateVector"]

        for item in state_vectors:
            if item["EPOCH"] == epoch:
                speed = instantaneous_speed(item)
                if speed is not None:
                    return jsonify({"epoch": epoch, "speed": speed})
                else:
                    return jsonify({"error": "Failed to calculate speed"})
        return jsonify({"error": "Epoch not found"})
    except Exception as e:
        logging.error(f"Failed to fetch speed for epoch {epoch}: {e}")
        return jsonify({"error": "An error occurred"})

# Route: /epochs/<epoch>/location
@app.route('/epochs/<epoch>/location', methods=['GET'])
def get_epoch_location(epoch):
    try:
        iss_data = json.loads(rd.get("iss_data"))
        state_vectors = iss_data["ndm"]["oem"]["body"]["segment"]["data"]["stateVector"]

        for item in state_vectors:
            if item["EPOCH"] == epoch:
                lat, lon, alt = compute_location_astropy(item)
                if lat is not None and lon is not None and alt is not None:
                    geocoder = Nominatim(user_agent='iss_tracker')
                    geoloc = geocoder.reverse((lat, lon), zoom=15, language='en')
                    return jsonify({
                        "epoch": epoch,
                        "latitude": lat,
                        "longitude": lon,
                        "altitude": alt,
                        "geoposition": geoloc.address if geoloc else "Unknown"
                    })
                else:
                    return jsonify({"error": "Failed to compute location"})
        return jsonify({"error": "Epoch not found"})
    except Exception as e:
        logging.error(f"Failed to fetch location for epoch {epoch}: {e}")
        return jsonify({"error": "An error occurred"})

# Route: /now
@app.route('/now', methods=['GET'])
def get_now():
    try:
        iss_data = json.loads(rd.get("iss_data"))
        state_vectors = iss_data["ndm"]["oem"]["body"]["segment"]["data"]["stateVector"]

        closest_entry = find_closest_epoch(state_vectors)
        if closest_entry:
            speed = instantaneous_speed(closest_entry)
            lat, lon, alt = compute_location_astropy(closest_entry)
            if lat is not None and lon is not None and alt is not None:
                geocoder = Nominatim(user_agent='iss_tracker')
                geoloc = geocoder.reverse((lat, lon), zoom=15, language='en')
                return jsonify({
                    "epoch": closest_entry["EPOCH"],
                    "speed": speed,
                    "latitude": lat,
                    "longitude": lon,
                    "altitude": alt,
                    "geoposition": geoloc.address if geoloc else "Unknown"
                })
            else:
                return jsonify({"error": "Failed to compute location"})
        return jsonify({"error": "No data available"})
    except Exception as e:
        logging.error(f"Failed to fetch current ISS data: {e}")
        return jsonify({"error": "An error occurred"})

# Run the Flask app
if __name__ == '__main__':
    load_iss_data()
    app.run(host='0.0.0.0', port=5000, debug=True)
