import btk


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
        writer.SetInput(clone)
        writer.SetFilename(filename.split('.')[0] + '-K' + str(i+1) + '.c3d')
        writer.Update()

a=1