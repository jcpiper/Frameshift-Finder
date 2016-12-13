# Frameshift-Finder
A Django app to identify frame shift loci in bacteriophage genomes

# Python Version #
2.7.11

# Django Version #
1.10.3

# App Name #
fsfinder

For info on django's architechture please see this [django tutorial](https://docs.djangoproject.com/en/1.10/intro/tutorial01/).

To run the app on your local machine, you will need to create a django project to contain the app: (this repository currently only holds the app file tree, not an entire project). In your console, navigate to the directory you would like to work in and run `$ django-admin startproject [arbitrary project name]`. Next, build the app directory. Move into the project directory -- `$ cd [arbitrary project name]` and execute `$ python manage.py startapp fsfinder`. This will create the `fsfinder` directory within your django project. You can then fork this directory (click the `fork` button in the top right corner of this page) and clone the repository into your fsfinder directory. To do this, make sure you are in your fsfinder directory. Now navigate to your fork of this repository (there should be a link on the right hand side of your github homepage) and click the `clone or download` button in the top right corner. This will diplay a url, copy this url and paste it into the command `$ git clone [your clone link]` within your console (i.e. if we were to clone the original repo, it would look like `$ git clone https://github.com/jcpiper/Frameshift-Finder.git`). That's it! You should be good to go after executing these commands. 

Run app locally with `python manage.py runserver` from command line. Note, you must execute this command while in the django project directory (the parent directory of fsfinder). You cannot execute this command from the fsfinder directory.
