from ezc3d import c3d
from pyomeca import Markers, Analogs
import config as cfg
import importlib
import numpy as np

def normalize_emg(emg):
    """
    Funkcja do normalizacji sygnalu emg
    """
    emg_p = (
    emg.meca.band_pass(order=2, cutoff=[10, 425])
    .meca.center()
    .meca.abs()
    .meca.low_pass(order=4, cutoff=5, freq=emg.rate)
    .meca.normalize(ref=None, scale=1)
    .meca.time_normalize(n_frames=425)
    )
    return emg_p

def name_body(data_path, marker_list):
    """
    Funkcja dodajaca numer badanego do nazwy markerow (potrzebne do wczytywania konkretnych markerów z pliku)
    """
    s = data_path.split('\\')[-1]
    body = s.split('-')[3]+":"
    for i in range(len(marker_list)):
        marker_list[i] = body+marker_list[i]
    return marker_list

def data_markers(data_path, marker_list):
    """
    Funkcja wczytujaca markery z pliku i zwracajaca tablice z danymi
    """
    frame_rate = 425  #Liczba klatek uzyta do normalizacji 
    
    #Wczytanie, normalizacja
    for i in range(len(marker_list)):
        tmp_markers = Markers.from_c3d(data_path, usecols=[marker_list[i]])
        tmp_markers = tmp_markers.meca.time_normalize(n_frames=frame_rate)
        tmp_markers = tmp_markers.meca.to_wide_dataframe()
        data_markers = data_markers.join(tmp_markers)      
        
    cols = [c for c in data_markers.columns if c.lower()[:4] != 'ones']
    markers_dataframe = data_markers[cols]
    return markers_dataframe

def dataframe(data_path):
    mark_channels = plat_channels = ang_channels = forc_channels = mome_channels = []
    """
    Funkcja zwracajaca tablice o rozmiarze ( 425, 190) - zawiera wszystkie odczytane dane
    """
    mark_channels = name_body(data_path,cfg.marker_channels[:13])
    plat_channels = name_body(data_path,cfg.platform_channels[:])
    ang_channels = name_body(data_path,cfg.angles_channels[:])
    forc_channels = name_body(data_path,cfg.forces_channels[:])
    mome_channels = name_body(data_path,cfg.moments_channels[:])
    powe_channels = name_body(data_path,cfg.powers_channels[:])
 
    try:    
	# Wywolywanie po kolei funkcji do wczytwania dla danych w odpowiedniej kolejnosci, i laczenie ich w odpowiedniej kolejnosc 
        full_data = data_markers(data_path,mark_channels)
        full_data = full_data.join(data_markers(data_path,cfg.marker_channels[13:]))
        emg1 = Analogs.from_c3d(data_path, usecols=cfg.emg_channels[:9])
        emg2 = Analogs.from_c3d(data_path, usecols=cfg.emg_channels[9:])
        norm_emg1 = normalize_emg(emg1)
        norm_emg2 = normalize_emg(emg2)
        norm_emg1 = norm_emg1.meca.to_wide_dataframe()
        norm_emg2 = norm_emg2.meca.to_wide_dataframe()

        full_data = full_data.join(norm_emg1)
        full_data = full_data.join(norm_emg2)

        full_data = full_data.join(data_markers(data_path,plat_channels))
        full_data = full_data.join(data_markers(data_path,ang_channels))
        full_data = full_data.join(data_markers(data_path,forc_channels))    
        full_data = full_data.join(data_markers(data_path,mome_channels))

        data_powers_tmp = data_markers(data_path,powe_channels)
        cols = [c for c in data_powers_tmp.columns if c.lower()[:1] == 'z'] # Wybranie tylko kolumn ze wspolrzedna z (dla power wartosci sa tylko z osi z)
        data_powers = data_powers_tmp[cols]

        full_data = full_data.join(data_powers)
        return full_data
    except:
        print("Blad wystapil w pliku: ",data_path)


#Petla do wyszukania wszystkich plikow i wywolania dla nich funkcji tworzacej tablice danych
y = np.zeros((3789,425,190)) 
for path in Path(r'Z:\baza').rglob('*.c3d'):
#     print("Sciezka: ",path, i)
    n = dataframe(str(path))
    y[i] = n # dodanie wczytanych danych do tablicy wynikowej
   
