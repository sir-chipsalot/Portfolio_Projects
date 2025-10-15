# ISS Tracker Application

## Description
This project is a Flask-based web application that tracks the International Space Station (ISS) using data provided by NASA. The application fetches the ISS's positional data, processes it, and provides various endpoints to query the ISS's location, speed, and other relevant information. The data is stored in Redis for quick access, and the application is containerized using Docker for easy deployment.

## Important Files
- `iss_tracker.py`: The main Flask application file that defines the routes and handles the logic for fetching and processing ISS data.
- `fetch_iss_tracker.py`: A script that fetches ISS data from NASA's API, parses it, and stores it in Redis and a JSON file.
- `Dockerfile`: Defines the Docker image for the Flask application.
- `docker-compose.yml`: Orchestrates the deployment of the Flask application and Redis database using Docker Compose.
- `Iss_data.json`: A JSON file that stores the fetched ISS data locally.
- `requirements.txt`: Lists the Python dependencies required for the application.

## Data Source
The ISS data is sourced from NASA's public API:
- [NASA ISS Trajectory Data](https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml)

## Deployment with Docker Compose
To deploy the application using Docker Compose, follow these steps:

1. **Clone the repository**:
   After cloning the repository run this command
    ```bash
    cd ~
    cd iss-tracker
    ```
2. **Remove data directory and then create it again**
    ```bash
    rm -r data
    mkdir data
    ```
3. **RUN this docker command IN THE ISS directory**
    ```bash
    docker compose up -d
   ```
4. **Enjoy your Iss tracker**
    ```bash
    you can run these commands:
    curl localhost:5000/epochs: Return entire data set including epochs, locations, and speeds
    curl localhost:5000/epochs?limit=int&offset=int: Returns a certain amount of epochs starting from the offset
    curl localhost:5000/epochs/<epoch>: Returns the data of a specific epoch set
    curl localhost:5000/epochs/<epoch>/speed: Returns the speed of a specific epoch
    curl localhost:5000/epochs/<epoch>/location: Returns the latitude, longitude, altitude, and geoposition for a specic epoch
    curl localhost:5000/now: Return instantaneous speed, latitude, longitude, altitude, and geoposition for the Epoch that is nearest in time
```
5. **Close the dockerfile**
```
    docker compose down
```
