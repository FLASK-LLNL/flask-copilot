FROM python:3.11-bookworm

RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
RUN . $HOME/.nvm/nvm.sh && nvm install 22

WORKDIR /app
COPY requirements.txt /app
COPY flask-app/. /app
COPY charge-backend /app/charge-backend

RUN . $HOME/.nvm/nvm.sh && \
    npm install && \
    npm run build

RUN pip install -r requirements.txt
RUN git clone --recursive https://github.com/FLASK-LLNL/ChARGe.git charge && \
    cd /app/charge && \
    pip install . && \
    charge-install --extras aizynthfinder --extras autogen --extras rdkit --extras chemprice

COPY mock_server.py /app
COPY dockerscripts/launch_servers.sh /app

RUN chmod -R g+rx /app

ENV FLASK_APPDIR=/app/dist

ARG AZF_PATH=./aizynth
COPY ${AZF_PATH}/. /aizynth

EXPOSE 8001
CMD ["--host", "0.0.0.0", "--port", "8001", "--workers", "8", "charge_server:app", "--json_file", "known_molecules.json", "--config", "/aizynth/config.yml", "--backend", "openai", "--model", "gpt-5-nano"]
ENTRYPOINT ["/app/launch_servers.sh"]
