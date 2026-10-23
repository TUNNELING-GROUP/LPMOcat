"""
URL configuration for the LPMOcat web application (lpmocat_web/).

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.shortcuts import redirect
from django.templatetags.static import static
from django.urls import path, include
from django.views.generic import RedirectView
from catalogue import views

urlpatterns = [
    path('', views.home, name='home'),
    # How to use and Acknowledgement are sections of the About page now; the old addresses keep working
    path('howtouse/', RedirectView.as_view(url='/#how-to-use', permanent=True)),
    path('acknowledgement/', RedirectView.as_view(url='/#acknowledgement', permanent=True)),
    path('catalogue/', include('catalogue.urls')),
    # path('catalogue/', views.page_in_progress, name='catalogue_disabled'), #update catalogue - test page
    path('admin/', admin.site.urls),
    # browsers and crawlers that ask for /favicon.ico without reading the pages' <link rel="icon">
    path('favicon.ico', lambda request: redirect(static('images/favicon.png'), permanent=True)),
    path('django_plotly_dash/', include('django_plotly_dash.urls'))
]
