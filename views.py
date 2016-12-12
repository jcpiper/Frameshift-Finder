from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.shortcuts import render
from .forms import UploadFileForm

# Create your views here.

# Algorithm will go here

def index(request):
	# render home page
	return render(request, 'fsfinder/index.html')
	
def upload(request):
	# handles file uploads
	if request.method == 'post':
		form = UploadFileForm(request.POST, request.FILES)
		if form.is_valid():
			process_file(request.FILES['file'])
			return HttpResponseRedirect('fsfinder/results.html')
		else:
			form = UploadFileForm()
		return render(request, 'fsfinder/index.html', {'form':form})
		
def process_file(file):
	# writes upload to temporary file
	with open('fsfinder/.fasta', 'w+') as destination:
		for chunk in file.chunks():
			destination.write(chunk)

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