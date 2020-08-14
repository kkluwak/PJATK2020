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

def wykresy_markerow(path):
    c = c3d(path)
    n_markers = ["LSHO","LELB","LWRA","RSHO","RELB","RWRA","RASI","RKNE","RANK"] # list waznych markerow
    axes = ["x","y","z"]
    body = path.split('-')[3]+":"
    p,k = dp.read_labels(path,200)
    for mark in markers:
        n = c['parameters']['POINT']['LABELS']['value'][0:44].index(body+mark) 
        for i in range(3):
            for j in range(len(p)):
                plt.plot(c['data']['points'][i][n][p[j]:k[j]])
            plt.title(axes[i])
            plt.show()
            
def nowy_czas(numer_markera,ev,markers):
    s=np.zeros(len(ev[0]))
    k=np.zeros(len(ev[0]))

    for i in range(len(ev[0])):
        output_difference=np.diff(markers[1][numer_markera][ev[0][i]:ev[1][i]])
        plt.plot(output_difference)

        dz=max(output_difference)*0.2
        s[i]=np.argmax(output_difference>dz)-40
        k[i]=len(output_difference) - np.argmax(output_difference[::-1]>dz)+40

        if s[i]<0:
            s[i]=0
        if k[i]>ev[1][i]:
            k[i]=ev[1][i]
        print('s',s[i],'k',k[i])
    return [s,k]

def przesuwanie_wykresow(ev,os,numer_markera,s,k,markers):
    
    evi=np.zeros((len(ev),len(ev[0])))
    for i in range(len(ev[0])):
        k.astype(int)
        s.astype(int)
        evi[1][i]=ev[0][i]+k[i]
        evi[0][i]=ev[0][i]+s[i]

    for i in range(len(evi[0])):
        print('start',s[i],'koniec',k[i])
        markers[os][numer_markera][int(evi[0][i]):int(evi[1][i])]=(markers[os][numer_markera][int(evi[0][i]):int(evi[1][i])]-min(markers[os][numer_markera][int(evi[0][i]):int(evi[1][i])]))/(max(markers[os][numer_markera][int(evi[0][i]):int(evi[1][i])])-min(markers[os][numer_markera][int(evi[0][i]):int(evi[1][i])]))

        t_konc=100
        dl_ciagu=int(evi[1][i])-int(evi[0][i])
        x=np.linspace(0,t_konc, dl_ciagu)
        plt.plot(x, markers[os][numer_markera][int(evi[0][i]):int(evi[1][i])])