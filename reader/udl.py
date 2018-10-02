# -*- coding: utf-8 -*-
#Unifed Data Loader

from numpy import array, double
import gzip
import bz2

delimlist = ['', ' ', ',', '\t', ';']

def fread(filename):
	'''This function read data from abstract file'''
	if filename.lower().endswith('.gz'):
		fp = gzip.open(filename, 'r')
	elif filename.lower().endswith('.bz2'):
		fp = bz2.BZ2File(filename, 'r')
	else:
		fp = open(filename, 'r')

	#ищем начало данных
	isdeccomma = False
	cnt = 0
	cntlim = 100
	inddelim = 0

	fstr = fp.readline()
	while cnt < cntlim and fstr:
		inddelim = find_delimiter( fstr)
		if inddelim < 0:
			inddelim = find_delimiter( fstr.replace(',','.'))
			if inddelim < 0:
				fstr = fp.readline()
				cnt += 1
				continue
			else:
				isdeccomma = True
				break
		else:
			break

	print("isdeccomma", isdeccomma)
	#если заголовок более 100 строк, что-то не так
	if cnt == cntlim or not fstr:
		return None

	#данные опознаны, считаем колонки
	if not isdeccomma:
		colcnt = check_columns4float( fstr, delimlist[inddelim])
	else:
		colcnt = check_columns4float( fstr.replace(',','.'), delimlist[inddelim])
	
	#если колонки не найдены, что-то не так
	if colcnt == 0:
		return None
	print("colcnt", colcnt)
	print("delimeter", delimlist[inddelim])

	d = []
	for i in xrange(colcnt):
		d.append(array([], dtype = double))

	cnt = 0
	#if isdeccomma:
	#	fstr = fstr.replace(',','.')
	#lstr = fstr.strip().split(delimlist[inddelim])
	#lstr = rm_empty(lstr)
	#ds = array([lstr], dtype='|S')
	filedata = fp.readlines();
	#filedata.insert(0, fstr)

	'''
	for fstr in filedata:
		if fstr.strip() == '': 
			break
		if isdeccomma:
			fstr = fstr.replace(',','.')
		lstr = fstr.strip().split(delimlist[inddelim])
		lstr = rm_empty(lstr)

		ds = append(ds, array(lstr, dtype='|S'))

		#for i in xrange(len(lstr)):
		#	d[i] = append( d[i], float(lstr[i]))

		fstr = fp.readline()
		cnt += 1
	'''
	if isdeccomma:
		filedata = [s.replace(',', '.') for s in filedata]

	#strip by column num if needed
	if delimlist[inddelim] is not '':
		if len(filedata[0].split(delimlist[inddelim])) != colcnt:
			filedata = [ delimlist[inddelim].join( s.split(delimlist[inddelim])[:colcnt]) for s in filedata]
		filedata = array( delimlist[inddelim].join( filedata ).split(delimlist[inddelim]))
	else:
		filedata = array( filedata)

	fp.close()

	#strip last empty strings if needed
	while filedata[-1].strip() is '':
		filedata = filedata[:-1]

	d = filedata.astype(double).reshape(-1, colcnt).T

	return d

def rm_empty(l):
	while 1:
		try:
			l.remove('')
		except ValueError:
			break
	return l

def check_columns4float( dstr, dlt):
	if dlt is not '':
		sl = dstr.split(dlt)
	else:
		sl = [dstr]
	cnt = 0;
	#remove all empty elements
	sl = rm_empty(sl)
	#проверка на преобразование
	for s in sl:
		try:
			float(s)
			cnt += 1
		except ValueError:
			break
	return cnt

def check_delimiter4float( dstr, dlt):
	if dlt is not '':
		sl = dstr.split(dlt)
	else:
		sl = [dstr]
	cnt = 0;
	#remove all empty elements
	while 1:
		try:
			sl.remove('')
		except ValueError:
			break
	if len(sl) > 1:
		for s in sl:
			try:
				float(s)
				cnt += 1
			except ValueError:
				if s == sl[-1]: #ошибка в последнем столбике допустима
					cnt+=1
				else:
					return False
		if cnt == len(sl):
			return True
	else:
		try:
			float(sl[0])
			return True
		except ValueError:
			return False

def find_delimiter( s):
	for i in xrange(len(delimlist)):
		if check_delimiter4float( s, delimlist[i]):
			return i
	return -1

