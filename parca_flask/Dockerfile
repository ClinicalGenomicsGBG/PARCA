#https://www.digitalocean.com/community/tutorials/how-to-build-and-deploy-a-flask-application-using-docker-on-ubuntu-18-04

FROM tiangolo/uwsgi-nginx-flask:python3.8
MAINTAINER Pernilla Ericsson
ENV STATIC_URL /templates
ENV STATIC_PATH /var/www/app/templates
COPY ./requirements.txt /var/www/requirements.txt
RUN pip install -r /var/www/requirements.txt
