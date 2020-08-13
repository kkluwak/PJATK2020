import btk
from ezc3d import c3d
import numpy as np
import matplotlib.pyplot as plt

def cropp_c3dfile(eventsFrame, filename):
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
        writer.SetFilename(filename.split('.')[0] + '-K' + str(i+1) + '.c3d')
        writer.Update()

def read_labels(data_path,frame_rate):
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

def rysowanie_marker(data,point,p,k):
    for j in range(3):
        for i in range(10):
            plt.plot(data[j][point][p[i]:k[i]])
        plt.show()
    
a=1500
    
