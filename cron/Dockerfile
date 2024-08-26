FROM python:3.9-slim

RUN apt-get update && apt-get -y install cron

WORKDIR /root/
COPY cron/requirements.txt requirements.txt
COPY cron/script.sh script.sh

RUN mkdir /root/logs

# Give execution rights on the cron job
RUN chmod 755 script.sh

RUN pip install -Ur requirements.txt

COPY src/ src/
COPY resources/ resources/

# Copy hello-cron file to the cron.d directory
COPY cron/get_data-cron /etc/cron.d/get_data-cron

# Give execution rights on the cron job
RUN chmod 0644 /etc/cron.d/get_data-cron

# Apply cron job
RUN crontab /etc/cron.d/get_data-cron

# Create the log file to be able to run tail
RUN touch /var/log/cron.log

# Run the command on container startup
CMD ["cron","&&","tail","-f","/var/log/cron.log"]
