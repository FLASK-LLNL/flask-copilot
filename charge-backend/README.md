## Running the ChARGe Backend

1. **Install Dependencies**
    - Check the `requirements.txt` file for all necessary dependencies and install them using pip:
      ```bash
      pip install -r requirements.txt
      ```
2. **Start the MCP Server**
    - Running the `Lead Molecule Optimization` tools server:
      ```bash
      python lmo_tool_server.py
      ```
    - Run the `Reaction Prediction + FlaskV2` server (coming soon):
      ```bash
      python flaskv2_server.py --port 6000
      ```
    
3. **Run the Charge Backend**
    - Start the ChARGe backend server:
      ```bash
      python charge_server.py  --json_file known_molecules.json --server-urls <LMO_TOOL_SERVER_URL> <FLASKV2_SERVER_URL> 
      ```
      