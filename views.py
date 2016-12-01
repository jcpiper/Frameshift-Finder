from django.shortcuts import render

# Create your views here.

# Algorithm will go here

def index(request):
	# render home page
	return render(request, 'fsfinder/index.html')

def file_check(request):
	# check input file for FASTA format
	# if correct format, pass sequence(s) onto calculate view
	# if not, render an error page/allow user to submit another file
	valid = True # set to false if any format errors are found
	if valid:
		calculate(request)
	else:
		context = {
			'error': True,
			# could add in a list of errors here
		}
		return render(request, 'fsfinder/index.html', context)

def calculate(request):
	# Algorithm goes here
	
	# Another file to build DFA?
	return render(request, 'fsfinder/results.html')