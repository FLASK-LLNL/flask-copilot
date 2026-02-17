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
   - Running the `Retrosynthesis` tools server:
     ```bash
     python retro_tool_server.py --port <PORT> --config <CONFIG_FILE>
     ```
     - The `--config` argument should point to the configuration `yaml` file for `AiZynthFinder`.
   - Run the `Reaction Prediction + FlaskV2` server:
     ```bash
     python flaskv2_server.py --port 6000
     ```

3. **Run the Charge Backend**
   - Start the ChARGe backend server:
     ```bash
     python charge_server.py  --json_file known_molecules.json --server-urls <LMO_TOOL_SERVER_URL>/sse <FLASKV2_SERVER_URL>/sse <RETRO_TOOL_SERVER_URL>/sse
     --config_file <CONFIG_FILE>
     ```
     - The `--json_file` argument should point to a JSON file containing known molecules.
     - The `--server-urls` argument should include the URLs of the running MCP servers
     - The `--config_file` argument should point to the configuration `yaml` file for `AiZynthFinder`.
