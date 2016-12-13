from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.shortcuts import render
from .forms import UploadFileForm
import os
import subprocess
# import os.popen*

# def index(request):
	# render home page
	# return render(request, 'fsfinder/index.html')
	
def upload(request):
	# handles file uploads
	if request.method == 'POST':
		form = UploadFileForm(request.POST, request.FILES)
		
		if form.is_valid():
			process_file(request.FILES['file'])
			return HttpResponseRedirect('results.html')
		else:
			form = UploadFileForm()
			return render(request, 'fsfinder/index.html', {'form':form})
	else:
		form = UploadFileForm()
		return render(request, 'fsfinder/index.html', {'form':form})

def process_file(file):
	# writes upload to temporary file
	with open('fsfinder/upload.fasta', 'wb+') as destination:
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
	file = open('fsfinder/upload.fasta', 'r')
	if file.readline() == '':
		filestatus = 'No file found'
	else:
		filestatus = subprocess.Popen(["python", "fsfinder/FrameshiftFinder.py"], stdout=subprocess.PIPE).communicate()[0] # a little python black magic to access output of our algorithm
		output = []
		output.append(filestatus[1:38])
		cleanOutput = filestatus[38:]
		while len(cleanOutput) > 36:
			output.append(cleanOutput[5:42])
			cleanOutput = cleanOutput[42:]
			
	return render(request, 'fsfinder/results.html', {'filestatus':filestatus, 'candidates':output})