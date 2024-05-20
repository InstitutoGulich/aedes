FROM ubuntu:20.04
ENV PYTHONUNBUFFERED 1

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends pip gcc g++ \
    python3 python3-dev git apache2 libapache2-mod-python python3-venv

#RUN useradd --create-home appuser
#WORKDIR /home/appuser
#USER appuser

WORKDIR /var/www
COPY web/wrapper.sh wrapper.sh

WORKDIR /root/
COPY web/requirements.txt requirements.txt

RUN python3.8 -m venv venv
ENV PATH="/root/venv/bin:$PATH"

RUN pip install -Ur requirements.txt

RUN git clone https://github.com/pybind/pybind11.git
RUN cp pybind11/include/pybind11 /usr/include/  -r

RUN git clone https://gitlab.com/libeigen/eigen.git
RUN cp eigen/Eigen /usr/include/  -r

COPY resources/ resources/
COPY src/ src/
RUN g++ -std=c++17 -Wall -O3 -march=native -shared -fPIC -I/usr/include/python3.8 src/cpp/otero_precipitation_wrapper.cpp -o src/otero_precipitation_wrapper.so

RUN a2enmod cgid
RUN service apache2 restart
COPY web/cgi-enabled.conf /etc/apache2/conf-available/cgi-enabled.conf

RUN ln -s /root/src/app /var/www/html/app
RUN chmod -R +x /root/
RUN a2enconf cgi-enabled
RUN sed -i.bak 's/#AddDefaultCharset\ UTF-8/AddDefaultCharset\ UTF-8/' /etc/apache2/conf-available/charset.conf
RUN chmod 705 src/app/cgi-bin/iface.py
RUN service apache2 restart


RUN mkdir src/app/outputs
RUN chmod 777 src/app/outputs

EXPOSE 80

CMD /usr/sbin/apache2ctl -D FOREGROUND

#docker build -t dengue .
#docker run -d -p 8989:80 -v $(pwd)/data:/data --name=dengue_container dengue
