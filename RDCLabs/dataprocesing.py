import btk
from ezc3d import c3d
import numpy as np
import matplotlib.pyplot as plt

def cropp_c3dfile(eventsFrame, filename, destiny):
    """
    Funkcja oddzielajaca pojedyncze ruchy w odrebne pliki na podstawie danych o markerach
    
    Przy wyowolaniu nalezy podac poczatek i koniec wycinka w formacie [[a,b],[a,b],...],
    sciezke pliku, ktory zostanie pociety oraz sciezke, do ktorej zostana zapisane wyodrebnione czesci
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
    Funkcja zwraca tablice [p, d], w której są zapisane czasy eventow oznaczających
    przyjecie postawy poczatkowej.
   
    Wywolanie funkcji wymaga podania dwoch parametrow, pierwszy - sciezka do pliku c3d, drugi - frame rate
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

#Funkcja do obliczania nowego punktu startowego (s) i końcowego (k). Podaje się na wejściu numer markera według którego to liczymy (sugerowana prawa dłoń lub prawa stopa)
#ev to początek i koniec ruchu na podstawie eventów (to co zwraca funkcja read_labels), markers to współrzędne markerów (Markers.from_c3d(data_path, prefix_delimiter=":"))
def nowy_czas(numer_markera,ev,markers):
    s=np.zeros(len(ev[0]))
    k=np.zeros(len(ev[0]))
#liczymy rozniczkę dla każdego uderzenia (od eventu postawy początkowej do eventu postawy początkowej
    for i in range(len(ev[0])):
        output_difference=np.diff(markers[1][numer_markera][ev[0][i]:ev[1][i]])
        #plt.plot(output_difference)

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
        #print('s',s[i],'k',k[i])
        # 
    return [s,k]

#Funkcja do przesuwania przyciętych wykresów tak, żeby można je było na siebie nałożyć
def przesuwanie_wykresow(ev,os,numer_markera,s,k,markers):
    
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
            #print('start',s[i],'koniec',k[i])

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
    c = c3d(path)
    n_markers = ["LSHO","LELB","LWRA","RSHO","RELB","RWRA","RASI","RKNE","RANK"] # list waznych markerow
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