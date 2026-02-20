# FLASK Copilot Web UI

This is a Web UI for the FLASK Copilot, which presents computed molecules and
their properties. It can be used for reaction prediction, lead molecule
optimization, and other custom prompts.

The FLASK Copilot consists of a React application as a frontend, and a Python
WebSocket-powered server as the backend.

## Installing

- Backend:

  - Make sure you have Python installed with a virtualenv.
  - Go to the main folder and run `pip install -r requirements.txt`
  - Alternativey: To install the package, clone the repository and run:

  ```bash
  pip install -e .[all]
  flask-copilot-install --extras all
  ```

- Frontend:
  - Install `npm`
  - `cd` into the `flask-app` folder and run `npm install`
  - Go into the `flask-app` directory and run `npm start dev` for development
    work or `npm run build` for a production build of the app.

## Running

To run FLASK Copilot, both the frontend and the backend need to run. The backend
will also serve the frontend web UI on the same port, if `npm run build` was
run. If this is not the case (e.g., with `npm start dev`), the backend still
needs to run. A server that creates mock data will run with
`python mock_server.py`.

Note: if the server was not running when the web UI started, click the blinking
red dot on the top right side to reconnect.

## License

Copyright (c) 2025, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

SPDX-License-Identifier: Apache-2.0

LLNL-CODE-2006345
