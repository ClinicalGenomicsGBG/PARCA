# Maintainer Pernilla Ericsson
from flask import Flask
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.wrappers import Response

app = Flask(__name__)
app.wsgi_app = DispatcherMiddleware(
    Response('Not Found', status=404),
    {'/parca': app.wsgi_app}
)
from app import views
