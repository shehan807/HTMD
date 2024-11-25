import sys
import os
import yaml

# Function to load the config file and set environment variables
def load_config(config_file):
    try:
        # Read the YAML config file
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)

        # Set each item in the config as an environment variable
        for key, value in config.items():
            os.environ[key] = str(value)  # Convert everything to string for environment variable compatibility

    except FileNotFoundError:
        print(f"Error: Config file '{config_file}' not found.")
        sys.exit(1)
    except yaml.YAMLError as exc:
        print(f"Error parsing the YAML config file: {exc}")
        sys.exit(1)

