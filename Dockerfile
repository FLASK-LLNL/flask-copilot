FROM python:3.11-bookworm

RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
RUN . $HOME/.nvm/nvm.sh && nvm install 22

ARG DATA_PATH=./data
COPY ${DATA_PATH}/. /data
RUN mkdir /data/db

WORKDIR /app
COPY requirements.txt /app
COPY flask-app /app/flask-app
COPY charge_backend /app/charge_backend
COPY externals /app/externals

# Build the components for the LC Conductor
WORKDIR /app/externals/lc_conductor/lcc_ui_components
RUN . $HOME/.nvm/nvm.sh && \
    npm install -g npm@latest && \
    npm install
    # && \
    # rm -rf node_modules/.vite dist .vite
RUN . $HOME/.nvm/nvm.sh && npm run build

# Switch back to the flask-app to build copilot
WORKDIR /app/flask-app
RUN . $HOME/.nvm/nvm.sh && \
    npm install -g npm@latest && \
    npm install && \
    rm -rf node_modules/.vite dist .vite
RUN . $HOME/.nvm/nvm.sh && npm run build

# Switch back to the top level to build everything else
WORKDIR /app/

RUN python -m venv /venv
RUN . /venv/bin/activate && \
    pip install -e externals/lc_conductor && \
    pip install -r requirements.txt && \
    git clone --recursive https://github.com/FLASK-LLNL/ChARGe.git charge.git && \
    cd /app/charge.git && \
    pip install . && \
    charge-install --extras aizynthfinder --extras autogen --extras rdkit --extras chemprice

COPY mock_server.py /app
COPY dockerscripts/launch_servers.sh /app

RUN chmod -R g+rwx /app /data

ENV FLASK_APPDIR=/app/flask-app/dist

ARG SERVER_VERSION=""
ENV SERVER_VERSION=${SERVER_VERSION}

EXPOSE 8001
ENTRYPOINT ["/app/launch_servers.sh"]
