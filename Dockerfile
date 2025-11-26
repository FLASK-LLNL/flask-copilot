FROM python:3.11-bookworm

RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
RUN . $HOME/.nvm/nvm.sh && nvm install 22

WORKDIR /app
COPY requirements.txt /app
COPY flask-app/. /app
COPY charge-backend /app/charge-backend

ARG AZF_PATH=./aizynth
COPY ${AZF_PATH}/. /aizynth

RUN . $HOME/.nvm/nvm.sh && \
    npm install -g npm@latest && \
    npm install && \
    npm run build

RUN python -m venv /venv
RUN . /venv/bin/activate && \
    pip install -r requirements.txt && \
    git clone --recursive https://github.com/FLASK-LLNL/ChARGe.git charge.git && \
    cd /app/charge.git && \
    pip install . && \
    charge-install --extras aizynthfinder --extras autogen --extras rdkit --extras chemprice

COPY mock_server.py /app
COPY dockerscripts/launch_servers.sh /app

RUN chmod -R g+rwx /app /aizynth

ENV FLASK_APPDIR=/app/dist

ARG SERVER_VERSION=""
ENV SERVER_VERSION=${SERVER_VERSION}

EXPOSE 8001
ENTRYPOINT ["/app/launch_servers.sh"]
