import requests
import xmltodict
import json
import logging
import redis

# Set up logging to show warnings and errors
logging.basicConfig(level=logging.DEBUG)

def fetch_iss_data():
    # Hardcode Redis connection details for local use
    redis_host = 'redis-db'
    redis_port = 6379

    # Try to get the ISS data from NASA's public API
    try:
        response = requests.get(url="https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml")
    except requests.RequestException as e:
        logging.error(f"Failed to retrieve ISS data: {e}")
        exit(1)  # Stop the script if data can't be fetched

    # Try to convert the XML response into a Python dictionary
    try:
        data = xmltodict.parse(response.text)
    except Exception as e:
        logging.error(f"Failed to parse XML data: {e}")
        exit(1)

    # Save the parsed data into a JSON file for later use
    final_data = json.dumps(data, indent=2)
    try:
        with open("Iss_data.json", "w") as file:
            file.write(final_data)
    except Exception as e:
        logging.error(f"Failed to write JSON data to file: {e}")
        exit(1)

    # Save the data to Redis
    try:
        rd = redis.Redis(host=redis_host, port=redis_port, db=0)
        rd.set("iss_data", final_data)
        logging.info("ISS data fetched, saved to Redis, and written to JSON file.")
    except redis.RedisError as e:
        logging.error(f"Failed to save ISS data to Redis: {e}")
        exit(1)

if __name__ == '__main__':
    fetch_iss_data()
