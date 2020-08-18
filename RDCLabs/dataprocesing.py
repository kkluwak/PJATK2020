import btk
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import importlib
from ezc3d import c3d
from pyomeca import Analogs
from pyomeca import Markers

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
		dz=max(output_difference)*0.2
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
        
def wykresy_markerow(path,markers):
	"""
    Funkcja do wyświetlania wykresów markerów. 
	
	Input:
	- path - ściezka do pliku c3d
	- markers - współrzędne markerów (Markers.from_c3d(path, prefix_delimiter=":"))
	
	Output:
	- Wykresy położenia markerów
	
	"""
	c = c3d(path)
	n_markers = ["LSHO","LELB","LWRA","RSHO","RELB","RWRA","RASI","RKNE","RANK"] # lista waznych markerow
	axes = ["x","y","z"]
	body = path.split('-')[3]+":"
	p,k = read_labels(path,200)
	for mark in markers:
		n = c['parameters']['POINT']['LABELS']['value'][0:44].index(body+mark) 
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
    - datapath - ścieżka dostępu do fpliku c3d
    
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
            plt.plot(emg_processed_event)
        subplot(1, 2, 2)
        plt.plot(emg_processed[num])
        plt.show()

def compare_events_average(folder_path, person, exer_num):
    """
    Funkcja wyświetlająca uśrednioną prace mięsni dla danego świczenia i aktora.
    
    Input:
    - folder_path - ścieżka dostępu do folderu z wszystkimi nagraniami
    - person - nazwa aktora do wczytania
    - exer_num - numer ćwiczenia do wczytania
    
    Output:
    - Wykresy średnich przebiegów dla danego ćwiczenia
    
    """
    
    muscles_names2 = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]
    cons1="\*\*-E0"
    cons2="-*.c3d"
    path=folder_path+person+cons1+exer_num+cons2
     
    aver_arr_all=np.zeros((16,1000))     
    
    for file in glob.glob(path,recursive = True):
        print(file)
        emg_processed=emg_full_preproces(file)

        aver_arr=np.zeros((16,1000))  
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
                    aver_arr[num][t]=aver_arr[num][t]+time_normalized.values[t]

            aver_arr[num]=aver_arr[num]/10
            time=np.linspace(1,1000,1000)

            for t2 in range(1000):
                aver_arr_all[file_num][t2]=aver_arr_all[file_num][t2]+time_normalized.values[t2]
            file_num=file_num+1;        
    
    for num in range(16):
        subplot(1, 1, 1)
        plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=2, 
                    top=0.7, 
                    wspace=0.25, 
                    hspace=0.35)
        aver_arr_all[num]=aver_arr_all[num]/5
        plt.plot(time,aver_arr_all[num])     
        plt.title(muscles_names2[num])
        plt.show()
    print(aver_arr_all)	

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
	
    muscles_names2 = ["Czworoboczny grzbietu L","Trójgłowy ramienia L", "Dwugłowy ramienia L", "Prostownik nadgarstka L","Skośny brzucha L", "Pośladkowy średni L","Czworogłowy uda L", "Brzuchaty łydki L","Czworoboczny grzbietu P","Trójgłowy ramienia P", "Dwugłowy ramienia P", "Prostownik nadgarstka P","Skośny brzucha P", "Pośladkowy średni P","Czworogłowy uda P", "Brzuchaty łydki P"]
    cons1="\*\*-E0"
    cons2="-*.c3d"
    path=folder_path+person+cons1+exer_num+cons2
     
    aver_arr_all=np.zeros((16,1000))     
    
    for file in glob.glob(path,recursive = True):
        print(file)
        emg_processed=emg_full_preproces(file)
        
                
        aver_arr=np.zeros((16,1000))  
        file_num=0

        p,d=read_labels(file, 1000)
        ev=[p,d]
        for num in range(16):
            s,k=nowy_czas_analog(p,d,emg_processed[num])
          
            for i in range(len(p)):               
                emg_processed_event=emg_processed[num][(p[i]+s[i].astype(int)):(d[i]+k[i].astype(int))]
                emg_processed_event2 = (
                emg_processed_event.meca.normalize(scale=1)                
         )                                           
                time_normalized=emg_processed_event2.meca.time_normalize(n_frames=1000)                
    
                for t in range(1000):
                    aver_arr[num][t]=aver_arr[num][t]+time_normalized.values[t]

            aver_arr[num]=aver_arr[num]/10
            time=np.linspace(1,1000,1000)

            for t2 in range(1000):
                aver_arr_all[file_num][t2]=aver_arr_all[file_num][t2]+time_normalized.values[t2]
            file_num=file_num+1;
            
    for num in range(16):
        subplot(1, 1, 1)
        plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=2, 
                    top=0.7, 
                    wspace=0.25, 
                    hspace=0.35)
        aver_arr_all[num]=aver_arr_all[num]/5
        plt.plot(time,aver_arr_all[num])     
        plt.title(muscles_names2[num])
        plt.show()
