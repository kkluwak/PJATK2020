import btk
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import importlib
from ezc3d import c3d
from pyomeca import Markers, Analogs

from matplotlib.pyplot import subplot

def cropp_c3dfile(eventsFrame, filename, destiny):
	"""
	Funkcja oddzielajaca pojedyncze ruchy w odrebne pliki na podstawie danych o markerach.
	
	Input:
	-eventsFrame - poczatek i koniec wycinka w formacie [[a,b],[a,b],...]
	-filename - sciezka pliku do podzielenia
	-destiny - sciezka, do ktorej zostana zapisane wyodrebnione czesci
	
	Output:
	- Podzielone pliki c3d zawierajace dane o pojedynczym ruchu
	
	"""
	reader = btk.btkAcquisitionFileReader()
	reader.SetFilename(filename)
	reader.Update()
	acq = reader.GetOutput()
 
	writer = btk.btkAcquisitionFileWriter()
	
	for i in range(0, len(eventsFrame)):
		clone = acq.Clone();
		clone.ResizeFrameNumberFromEnd(acq.GetLastFrame() - eventsFrame[i][0] + 1)
		clone.ResizeFrameNumber(eventsFrame[i][1] - eventsFrame[i][0] + 1)
		clone.SetFirstFrame(eventsFrame[i][0])
		clone.ClearEvents()
		for e in btk.Iterate(acq.GetEvents()):
			if ((e.GetFrame() > clone.GetFirstFrame()) and (e.GetFrame() < clone.GetLastFrame())):
				  clone.AppendEvent(e)
		clone.SetFirstFrame(1)
		writer.SetInput(clone)
		writer.SetFilename(destiny + '\\' + (filename.split('\\')[-1]).split('.')[0]+ '-K' + str(i+1) + '.c3d')
		writer.Update()
	  
def read_labels(data_path,frame_rate):  
	"""
	Funkcja zwraca tablice [p, k], w której są zapisane czasy eventow oznaczających przyjecie postawy poczatkowej.
	
	Input:
	- data_path - sciezka do pliku c3d
	- frame_rate - częstotliwośc próbkowania danych w pliku

	Output:
	- [p,k] - tablice punktów startowych (s) i końcowych(k)

	"""
	
	c3d_to_compare= c3d(data_path)
	event = c3d_to_compare['parameters']['EVENT']['LABELS']['value']
	czas = np.around(c3d_to_compare['parameters']['EVENT']['TIMES']['value'][1]*frame_rate)
	eventy = [event, czas]
	
	eventy[0].index('Foot Strike')
	indxE = [i for i, x in enumerate(eventy[0]) if x == "Event"]
	indxFS = [i for i, x in enumerate(eventy[0]) if x == "Foot Strike"]

	CzasFS = np.zeros(len(indxFS))
	for i in range(len(indxFS)):
		CzasFS[i] = eventy[1][indxFS[i]]

	CzasE = np.zeros(len(indxE))
	for i in range(len(indxE)):
		CzasE[i] = eventy[1][indxE[i]]
	eventy[1].sort()

	p = []
	k = []
	for i in range(len(eventy[1])):
		if not i >= len(eventy[1])-2:
			pierwszy = eventy[1][i]
			drugi = eventy[1][i+1]
			trzeci = eventy[1][i+2]
			if pierwszy in CzasE:
				if drugi in CzasFS:
					if trzeci in CzasE:
						p.append(int(pierwszy))
						k.append(int(trzeci))
	return [p,k]

def nowy_czas_marker(numer_markera,ev,markers):
	"""
	Funkcja do obliczania nowego punktu startowego (s) i końcowego (k).
	
	Input:
	- numer_markera - numer markera według któreg liczymy nowe punkty (sugerowana prawa dłoń lub prawa stopa)
	- ev - początek i koniec ruchu na podstawie eventów (to co zwraca funkcja read_labels)
	- markers - współrzędne markerów (Markers.from_c3d(data_path, prefix_delimiter=":"))

	Output:
	- [s,k] - nowe punkty startowe (s) i końcowe(k)
	"""
	s=np.zeros(len(ev[0]))
	k=np.zeros(len(ev[0]))
	#liczymy rozniczkę dla każdego uderzenia (od eventu postawy początkowej do eventu postawy początkowej
	for i in range(len(ev[0])):
		output_difference=np.diff(markers[1][numer_markera][ev[0][i]:ev[1][i]])
		
		#ustalenie nowego startu i końca ruchu
		dz=max(output_difference)*0.2
		dx=min(output_difference)*0.9
		s[i]=np.argmax(output_difference>dz)-40
		k[i]=len(output_difference) - np.argmax(output_difference[::-1]>dx)+40

		#warunki, które mają zabezpieczać przed wyjściem za zakres pociętego nagrania
		if s[i]<0:
			s[i]=0
		if k[i]>ev[1][i]:
			k[i]=ev[1][i]
	return [s,k]
	
def nowy_czas_analog(p,d,analogs):
	"""
	Funkcja do obliczania nowego punktu startowego (s) i końcowego (k).
	
	Input:
	- p,d - początek i koniec ruchu na podstawie eventów (p,d - to co zwraca funkcja read_labels)
	- analogs - przegieg sygnału EMG dla wybranego mięśnia ((Analogs.from_c3d(datapath, usecols=muscles))[numer_mięśnia])

	Output:
	- [s,k] - nowe punkty startowe (s) i końcowe(k)
	
	"""
	
	ev=[p,d]
	s=np.zeros(len(ev[0]))
	k=np.zeros(len(ev[0]))
	#liczymy rozniczkę dla każdego uderzenia (od eventu postawy początkowej do eventu postawy początkowej
	for i in range(len(ev[0])):
		output_difference=np.diff(analogs[ev[0][i]:ev[1][i]])
		#ustalenie nowego startu i końca ruchu
		dz=max(output_difference)*0.3
		dx=min(output_difference)*0.9
		
		s[i]=np.argmax(output_difference>dz)
		k[i]=len(output_difference) - np.argmax(output_difference>dx)

		#warunki, które mają zabezpieczać przed wyjściem za zakres pociętego nagrania
		if s[i]<0:
			s[i]=0
		if k[i]>ev[1][i]:
			k[i]=ev[1][i] 
	return [s,k]

def przesuwanie_wykresow(ev,numer_markera,s,k,markers):
	"""
	Funkcja do wyświetlania wykresów markerów przesuniętych w fazie. 
	
	Input:
	- ev - początek i koniec ruchu na podstawie eventów (to co zwraca funkcja read_labels)
	- s,k - nowege punkty startowe i końcowe (zwraca je funkcja nowy_czas_marker)
	- numer_markera - określa dla którego markera chcemy wyświetlić wykres
	- markers - współrzędne markerów (Markers.from_c3d(data_path, prefix_delimiter=":"))
	
	Output:
	- Wykresy położenia zadanego markera
	
	"""
	#robimy z k i s inty, bo były z tym problemy
	#tworzymy nową tablicę zawierającą czasy startu i końca ruchu
	evi=np.zeros((len(ev),len(ev[0])))
	for i in range(len(ev[0])):
		k.astype(int)
		s.astype(int)
		evi[1][i]=ev[0][i]+k[i]
		evi[0][i]=ev[0][i]+s[i]
	
	#dla 3 osi robimy pętlę z robieniem wykresu
	for j in range(3):
		#dla ilości powtórzeń (zwykle 10) robimy pętlę żeby wyrzucało je na tym samym wykresie
		for i in range(len(evi[0])):

			#normalizacja
			markers[j][numer_markera][int(evi[0][i]):int(evi[1][i])]=(markers[j][numer_markera][int(evi[0][i]):int(evi[1][i])]-min(markers[j][numer_markera][int(evi[0][i]):int(evi[1][i])]))/(max(markers[j][numer_markera][int(evi[0][i]):int(evi[1][i])])-min(markers[j][numer_markera][int(evi[0][i]):int(evi[1][i])]))

			#operacja żeby rozciągnąć wykresy na tym samym okresie czasu (0-100)
			t_konc=100
			dl_ciagu=int(evi[1][i])-int(evi[0][i])
			x=np.linspace(0,t_konc, dl_ciagu)
			#plotowanie wykresu, w danej osi (ponieważ jest w pętli to zrobi się dla 3), dla danego numeru markera, od klatki startowej do końcowej
			plt.plot(x, markers[j][numer_markera][int(evi[0][i]):int(evi[1][i])])
		plt.show()
		
def wykresy_markerow(path,markers = ["LSHO","LELB","LWRA","RSHO","RELB","RWRA","RASI","RKNE","RANK"]):
	"""
	Funkcja do wyświetlania wykresów markerów. 
	
	Input:
	- path - ściezka do pliku c3d
	- markers - lista markerow do wyswietlenia
	
	Output:
	- Wykresy położenia markerów
	
	"""
	c = c3d(path)
	# markers = ["LSHO","LELB","LWRA","RSHO","RELB","RWRA","RASI","RKNE","RANK"] # lista waznych markerow
	axes = ["x","y","z"]
	body = path.split('-')[3]+":"
	p,k = read_labels(path,200)
	for mark in markers:
		n = c['parameters']['POINT']['LABELS']['value'].index(body+mark) 
		for i in range(3):
			for j in range(len(p)):
				plt.plot(c['data']['points'][i][n][p[j]:k[j]])
			plt.title(axes[i])
			plt.show()
			
def read_analog_allmuscles(datapath):
	"""
	Funkcja do wczytywania danych EMG dla wszystkich mięśni z pliku c3d. 
	
	Input:
	- datapath - ścieżka do pliku c3d z danymi
	
	Output:
	- emg - dane EMG z wczytanego pliku dla wszystkich mięśni
	
	"""
	muscles = ["Voltage.1","Voltage.2","Voltage.3","Voltage.4","Voltage.5","Voltage.6","Voltage.7","Voltage.8","Voltage.9","Voltage.10","Voltage.11","Voltage.12","Voltage.13","Voltage.14","Voltage.15","Voltage.16"]
	emg = Analogs.from_c3d(datapath, usecols=muscles)
	return emg
	
def rename_emg(emg):
	"""
	Funkcja do zmiany nazw kanałów EMG na nazwy odpowiadających im mięśni.
	
	Input:
	- emg - dane EMG z wczytanego pliku dla wszystkich mięśni z oryginalnymi nagłówkami (Analogs.from_c3d(datapath, usecols=muscles))
	
	Output:
	- emg - dane EMG z wczytanego pliku dla wszystkich mięśni ze zmienionymi nagłówkami
	
	"""
	muscles_names = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]
	emg['channel'] = muscles_names
	return emg
	
def show_emg_data(emg):
	"""
	Funkcja do wyświetlania danych EMG. 
	
	Input:
	- emg - dane EMG z pliku c3d (Analogs.from_c3d(datapath, usecols=muscles))
	
	Output:
	- Wyswietlenie wykresów emg 
	
	"""
	emg.plot(x="time", col="channel", col_wrap=3)
	
def normalize_emg(emg):
	"""
	Funkcja wykonująca normalizację danych EMG.
	
	Input:
	- emg - przegiegi sygnałów EMG (Analogs.from_c3d(datapath, usecols=muscles))

	Output:
	- Znormalizowane oraz przefiltrowane dane EMG

	"""
	emg_p = (
	emg.meca.band_pass(order=2, cutoff=[10, 425])
	.meca.center()
	.meca.abs()
	.meca.low_pass(order=4, cutoff=5, freq=emg.rate)
	.meca.normalize(ref=None, scale=1)
	)
	return emg_p 

def emg_full_preproces(datapath):
	"""
	Funkcja wczytująca dane EMG z pliku c3d oraz wstępnie je przetwarzająca.
	
	Input:
	- datapath - ścieżka dostępu do pliku c3d
	
	Output:
	- normalised_emg - znormalizowane dane EMG
	
	"""

	emg_data = read_analog_allmuscles(datapath)
	normalised_emg=normalize_emg(emg_data)
	return normalised_emg
	
def show_events(data_path):
	"""
	Funkcja wyświetla wykresy: po lewej ruchy nałozone na siebie w czasie, po prawej odczyt pełnej scieżki napieć mięsni w czasie.
	
	Input:
	- data_path - ścieżka do pliku c3d z danymi EMG
   
	Output:
	- Wykresy nałożonych na siebie ruchów oraz pełnego przegiegu pracy mięsnia dla całego nagrania

	"""
	muscles_names = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]

	
	emg_processed = emg_full_preproces(data_path)
	p,d=read_labels(data_path, 1000) # p - moment rozpoczecia eventu, d - momen zakończenia eventu
	print(p,d)
	for num in range(16):
		subplot(1, 2, 1)
		plt.subplots_adjust(left=0.125,
					bottom=0.1, 
					right=2.8, 
					top=0.9, 
					wspace=0.25, 
					hspace=0.35)
		for i in range(len(p)):             
			emg_processed_event=emg_processed[num][p[i]:d[i]]
			plt.plot(emg_processed_event, label=i+1)
			plt.title(muscles_names[num]+" - pojedyncze ruchy")
			plt.legend(loc='upper right')
		subplot(1, 2, 2)
		plt.plot(emg_processed[num])
		plt.title(muscles_names[num]+" - pelen przebieg")
		
		plt.show()

def compare_events_average(folder_path, person, exer_num):
	"""
	Funkcja wyświetla wykresy: po lewej ruchy nałozone na siebie w czasie z przesunięciem w fazie, po prawej odczyt pełnej scieżki napieć mięsni w czasie.
	
	Input:
	- data_path - ścieżka do pliku c3d z danymi EMG
   
	Output:
	- Wykresy nałożonych na siebie ruchów z przesunięciem w fazie oraz pełnego przebiegu pracy mięśnia dla całego nagrania

	"""
	
	muscles_names2 = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]
	cons1="\*\*-E0"
	cons2="-*.c3d"
	path=folder_path+person+cons1+exer_num+cons2
	 
	aver_arr=np.zeros((6,16,1000))     
	time=np.linspace(1,1000,1000)
	for file in glob.glob(path,recursive = True):
		print(file)
		emg_processed=emg_full_preproces(file)
 
		file_num=0

		p,d=read_labels(file, 1000)
		
		for num in range(16):

			for i in range(len(p)):

				emg_processed_event=emg_processed[num][p[i]:d[i]]
				emg_processed_event2 = (
				emg_processed_event.meca.normalize(scale=1)                
		 )
				time_normalized=emg_processed_event2.meca.time_normalize(n_frames=1000)

				for t in range(1000):
					aver_arr[file_num][num][t]=aver_arr[file_num][num][t]+time_normalized.values[t]
		file_num=file_num+1
	aver_arr_all=np.zeros((16,1000)) 
	for plik in range(file_num): 
		for num in range(16):
			for t in range(1000):
				aver_arr_all[num][t]= aver_arr_all[num][t]+aver_arr[plik][num][t]
	
	for num in range(16):
			for t in range(1000):
				aver_arr_all[num][t]= aver_arr_all[num][t]/(10*(file_num+1))
				   
	
	for num in range(16):
		subplot(1, 1, 1)
		plt.subplots_adjust(left=0.125,
					bottom=0.1, 
					right=2, 
					top=0.7, 
					wspace=0.25, 
					hspace=0.35)
		#aver_arr_all[num]=aver_arr_all[num]/5
		plt.plot(time,aver_arr_all[num])     
		plt.title(muscles_names2[num])
		plt.show()
	#print(aver_arr_all)	

def compare_events_average_shifted(folder_path, person, exer_num):
	"""
	Funkcja wyświetlająca uśrednioną prace mięsni dla danego świczenia i aktora z przesunięciem ruchów w fazie.
	
	Input:
	- folder_path - ścieżka dostępu do folderu z wszystkimi nagraniami
	- person - Nazwa aktora do wczytania
	- exer_num - Nazwa ćwiczenia do wczytania
	
	Output:
	- Wykresy średnich przebiegów dla danego ćwiczenia z przesunięciem ruchów w fazie
	
	"""
	
	muscles_names = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]
	cons1="\*\*-E0"
	cons2="-*.c3d"
	path=folder_path+person+cons1+exer_num+cons2
	 
	aver_arr=np.zeros((6,16,1000))     
	time=np.linspace(1,1000,1000)
	file_num=0
	
	for file in glob.glob(path,recursive = True):
		if file.find("fail") == -1 and file.find("Cal") == -1:
			print(file)
			emg_processed=emg_full_preproces(file)

			p,d=read_labels(file, 1000)
			if (file_num==0):
				max_frame, frame_size=find_max_frame(p,d,emg_processed[9][p[0]:d[0]])
			ev=[p,d]
			for num in range(16):


				for i in range(len(p)):     
					s,k=find_new_start(p[i],d[i],emg_processed[num],max_frame,frame_size,i+1)
					emg_processed_event=emg_processed[num][(p[i]+s):(d[i]+k)]
					emg_processed_event2 = (
					emg_processed_event.meca.normalize(scale=1)                
			)                                           
					time_normalized=emg_processed_event2.meca.time_normalize(n_frames=2000)
					time_normalized=time_normalized[:1000].meca.time_normalize(n_frames=1000)
                    #time_normalized=emg_processed_event2.meca.time_normalize(n_frames=1000)

					for t in range(1000):
						aver_arr[file_num][num][t]=aver_arr[file_num][num][t]+time_normalized.values[t]
			file_num=file_num+1
	aver_arr_all=np.zeros((16,1000))    
	for plik in range(file_num): 
		for num in range(16):
			for t in range(1000):
				aver_arr_all[num][t]= aver_arr_all[num][t]+aver_arr[plik][num][t]
	
	for num in range(16):
			for t in range(1000):
				aver_arr_all[num][t]= aver_arr_all[num][t]/(10*(file_num+1))
	for num in range(16):
		subplot(1, 1, 1)
		plt.subplots_adjust(left=0.125,
					bottom=0.1, 
					right=2, 
					top=0.7, 
					wspace=0.25, 
					hspace=0.35)
		std_div=aver_arr_all[num].std()
		#aver_arr_all[num]=aver_arr_all[num]/5
		plt.axhline(std_div,color="gray",label='Standard deviation')
		plt.plot(time,aver_arr_all[num])     
		plt.title(muscles_names[num])
		plt.legend(loc='upper left')
		plt.show()
	return aver_arr_all
		
		
def events_average_shifted(path):
	"""
	Funkcja wyświetlająca uśrednioną prace mięsni dla pojedynczego nagrania.
	
	Input:
	- path - ścieżka dostępu do pliku
	
	Output:
	- 
	
	"""
	
	muscles_names = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]
	
	 
	aver_arr=np.zeros((16,1000))     
	time=np.linspace(1,1000,1000)
	

	emg_processed=emg_full_preproces(path)
	


	p,d=read_labels(path, 1000)
	max_frame, frame_size=find_max_frame(p,d,emg_processed[9][p[0]:d[0]])
	ev=[p,d]
	for num in range(16):
		
	  
		for i in range(len(p)):     
			s,k=find_new_start(p[i],d[i],emg_processed[num],max_frame,frame_size,i+1)
			emg_processed_event=emg_processed[num][(p[i]+s):(d[i]+k)]
			emg_processed_event2 = (
			emg_processed_event.meca.normalize(scale=1)                
	 )                                           
			time_normalized=emg_processed_event2.meca.time_normalize(n_frames=2000)
			time_normalized=time_normalized[:1000].meca.time_normalize(n_frames=1000)
			#time_normalized=emg_processed_event2.meca.time_normalize(n_frames=1000)
			
			for t in range(1000):
				aver_arr[num][t]=aver_arr[num][t]+time_normalized.values[t]

		
	aver_arr_all=np.zeros((16,1000))    
	for num in range(16):
		for t in range(1000):
			aver_arr_all[num][t]= aver_arr_all[num][t]+aver_arr[num][t]
	
	for num in range(16):
			for t in range(1000):
				aver_arr_all[num][t]= aver_arr_all[num][t]/(10)
	for num in range(16):
		subplot(1, 1, 1)
		plt.subplots_adjust(left=0.125,
					bottom=0.1, 
					right=2, 
					top=0.7, 
					wspace=0.25, 
					hspace=0.35)
		plt.plot(time,aver_arr_all[num])     
		plt.title(muscles_names[num])
		plt.show()

		
		
def show_events_norm_shifted(data_path):
	"""
	Funkcja wyświetlająca prace mięsni dla danego świczenia i aktora z przesunięciem ruchów w fazie.
	
	Input:
	- data_path - ścieżka dostępu do pliku c3d
	
	Output:
	- Wykresy przebiegów dla danego ćwiczenia z przesunięciem ruchów w fazie
	
	"""
	
	emg_processed = emg_full_preproces(data_path)

	muscles_names = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]

	p,d=read_labels(data_path, 1000)
	max_frame, frame_size=find_max_frame(p,d,emg_processed[9][p[0]:d[0]])
	for num in range(16):
		
		subplot(1, 2, 1)
		plt.subplots_adjust(left=0.125,
					bottom=0.1, 
					right=2.8, 
					top=0.9, 
					wspace=0.25, 
					hspace=0.35)

		for i in range(len(p)): 
					   
			
			s,k=find_new_start(p[i],d[i],emg_processed[num],max_frame,frame_size,i+1)
			
			
			emg_processed_event=emg_processed[num][(p[i]+s):(d[i]+k)]
			emg_processed_event2 = (
			emg_processed_event.meca.normalize(scale=1)                
			)
			time_normalized=emg_processed_event2.meca.time_normalize(n_frames=2000)
	#             if i==9:
	#                 time_normalized=time_normalized[:1000].meca.time_normalize(n_frames=(int)(1000*(len(emg_processed[9][p[9]:d[9]])/k)))
	#             else:
	#                 time_normalized=time_normalized[:1000].meca.time_normalize(n_frames=1000)
				
			time_normalized=time_normalized[:1000].meca.time_normalize(n_frames=1000)
				
			plt.plot(time_normalized, label=i+1)     
			plt.title(muscles_names[num])
			plt.legend(loc='upper left')

		subplot(1, 2, 2)
		plt.plot(emg_processed[num])
		plt.title(muscles_names[num])
		plt.show()
		
		
		
def find_max_frame(p,d,analogs):
	val_arr=[]
	for frame in range(len(analogs)):
		val_arr.append(analogs[frame].values)
	max_val=max(val_arr)
	frame_size=len(analogs)

	return [val_arr.index(max_val),frame_size]
	
	
def find_new_start(p,d,analogs,max_frame,frame_size,event_num):

	val_arr=[]
	analogs=analogs[p:d]
	analogs=analogs.meca.time_normalize(n_frames=frame_size)
	
	for frame in range(len(analogs)):
		val_arr.append(analogs[frame].values)
	max_val=max(val_arr)
	this_max_frame=val_arr.index(max_val)
	 
	s=p+(this_max_frame-p-max_frame)
	k=d+(this_max_frame-p-max_frame) 

	return [s,k]
	
def find_new_start_base(p,d,analogs,max_frame,frame_size,event_num):

	val_arr=[]
	analogs=analogs[p:d]
	
	for frame in range(len(analogs)):
		val_arr.append(analogs[frame].values)
	max_val=max(val_arr)
	this_max_frame=val_arr.index(max_val)
	 
	s=p+(this_max_frame-p-max_frame)
	k=d+(this_max_frame-p-max_frame) 

	return [s,k]




def average_shifted_filled(folder_path, person, exer_num):
	"""
	Funkcja wyświetlająca uśrednioną prace mięsni dla danego świczenia i aktora z przesunięciem ruchów w fazie.
	
	Input:
	- folder_path - ścieżka dostępu do folderu z wszystkimi nagraniami
	- person - Nazwa aktora do wczytania
	- exer_num - Nazwa ćwiczenia do wczytania
	
	Output:
	- Wykresy średnich przebiegów dla danego ćwiczenia z przesunięciem ruchów w fazie
	
	"""
	
	muscles_names = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]
	cons1="\*\*-E0"
	cons2="-*.c3d"
	path=folder_path+person+cons1+exer_num+cons2
	 
	aver_arr=np.zeros((6,16,1000))     
	time=np.linspace(1,1000,1000)
	file_num=0
	
	for file in glob.glob(path,recursive = True):
		print(file)
		emg_processed=emg_full_preproces(file)
		
			
		

		p,d=read_labels(file, 1000)
		if (file_num==0):
			max_frame, frame_size=find_max_frame(p,d,emg_processed[9][p[0]:d[0]])
		ev=[p,d]
		for num in range(16):
			
		  
			for i in range(len(p)):     
				s,k=find_new_start(p[i],d[i],emg_processed[num],max_frame,frame_size,i+1)
				emg_processed_event=emg_processed[num][(p[i]+s):(d[i]+k)]
				emg_processed_event2 = (
				emg_processed_event.meca.normalize(scale=1)                
		 )                                           
				time_normalized=emg_processed_event2.meca.time_normalize(n_frames=2000)
				time_normalized=time_normalized[:1000].meca.time_normalize(n_frames=1000)
				#time_normalized=emg_processed_event2.meca.time_normalize(n_frames=1000)
				
				for t in range(1000):
					aver_arr[file_num][num][t]=aver_arr[file_num][num][t]+time_normalized.values[t]
		file_num=file_num+1
		
	aver_arr_all=np.zeros((16,1000))    
	for plik in range(file_num): 
		for num in range(16):
			for t in range(1000):
				aver_arr_all[num][t]= aver_arr_all[num][t]+aver_arr[plik][num][t]
	
	for num in range(16):
			for t in range(1000):
				aver_arr_all[num][t]= aver_arr_all[num][t]/(10*(file_num+1))
	subplot(1, 1, 1)
	plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=5, 
                    top=0.25, 
                    wspace=0.25, 
                    hspace=0.35)
	#std_div=aver_arr_all[num].std()
	#aver_arr_all[num]=aver_arr_all[num]/5
	#plt.axhline(std_div,color="gray",label='Standard deviation')
	plt.plot(time,aver_arr_all[num])
	plt.fill_between(time,aver_arr_all[num])
	plt.title(muscles_names[num])
	#plt.legend(loc='upper left')
	plt.show()
        
        
def print_onset_offset(aver_arr_all,freq=1000, fun="mean", mult=1 ,above=None ,below=None ):  
    if above is None:
        above=freq/2
    if below is None:
        below=freq/2
    
    muscles_names = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]
    time=np.linspace(1,1000,1000)
    freq=1000
    for num in range(16):
        
        
        #onsets=detect_onset_new(aver_arr_all[num],threshold= aver_arr_all[num].mean(), n_above=freq / 4,n_below=freq / 4 )
        #print(onsets)
        if fun == "mean":
            onsets=detect_onset_new(aver_arr_all[num],threshold= aver_arr_all[num].mean()*mult, n_above=above,n_below=below )
        if fun == "std":
            onsets=detect_onset_new(aver_arr_all[num],threshold= aver_arr_all[num].std()*mult, n_above=above,n_below=below )

        
#         onsets =  aver_arr_all[num].meca.detect_onset(
#         threshold= aver_arr_all[num].mean(),   # mean of the signal 
            
#         #threshold2= emg[i].std(),
            
#         n_above=freq / 4,                     # we want at least 1/2 second above the threshold
#         n_below=freq / 4,                     # we accept point below threshold for 1/2 second
#         ) 
        subplot(1, 1, 1)
        for (start, end) in onsets:   
            plt.axvspan(start, end, color="b")           
         
        
        plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=5, 
                    top=0.25, 
                    wspace=0.25, 
                    hspace=0.35)
        aver_arr_all[num]=aver_arr_all[num]/5
        plt.plot(time,aver_arr_all[num], color="r")     
        plt.title(muscles_names[num]+" dla osoby B0446")
        plt.xlabel("Time [ms]")
        plt.ylabel("Signal [mV]")
        plt.show()
        
        
def detect_onset_new( x, threshold, n_above: int = 1, n_below: int = 0, threshold2: int = None, n_above2: int = 1):
    if x.ndim != 1:
        raise ValueError(
            f"detect_onset works only for one-dimensional vector. You have {x.ndim} dimensions."
        )
    if isinstance(threshold, xr.DataArray):
        threshold = threshold.item()
    if isinstance(threshold2, xr.DataArray):
        threshold2 = threshold2.item()

    x = np.atleast_1d(x.copy())
    x[np.isnan(x)] = -np.inf
    inds = np.nonzero(x >= threshold)[0]
    if inds.size:
        # initial and final indexes of almost continuous data
        inds = np.vstack(
            (
                inds[np.diff(np.hstack((-np.inf, inds))) > n_below + 1],
                inds[np.diff(np.hstack((inds, np.inf))) > n_below + 1],
            )
        ).T
        # indexes of almost continuous data longer than or equal to n_above
        inds = inds[inds[:, 1] - inds[:, 0] >= n_above - 1, :]
        # minimum amplitude of n_above2 values in x to detect
        if threshold2 is not None and inds.size:
            idel = np.ones(inds.shape[0], dtype=bool)
            for i in range(inds.shape[0]):
                if (
                    np.count_nonzero(x[inds[i, 0] : inds[i, 1] + 1] >= threshold2)
                    < n_above2
                ):
                    idel[i] = False
            inds = inds[idel, :]
    if not inds.size:
        inds = np.array([])
    return inds
        