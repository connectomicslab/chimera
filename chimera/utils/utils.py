import json
import os

# Loading the JSON file containing the available parcellations
def _load_parctype_json(json_file=None):
    """
    Load the JSON file containing the available parcellations

    Returns
    -------
    data : dict
        Dictionary containing the available parcellations
    """

    cwd = os.getcwd()
    serJSON = os.path.join(cwd, 'parcTypes.json')
    with open(serJSON) as f:
        data = json.load(f)

    return data