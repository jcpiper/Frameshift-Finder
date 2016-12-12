# urls
from django.conf.urls import url
from . import views

urlpatterns = [
	url(r'^$', views.upload, name='index'),
	url(r'^results', views.calculate, name='results'),
]