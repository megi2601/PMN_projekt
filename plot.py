import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

def params_dict(filename):
    d = dict()
    params = filename.split("_")[1:]
   # u = float(params[0][1:])
    V=float(params[0][2:])
    N=int(params[1][1:])
    dt = float(params[2][2:])
    s = int(params[3][1:])
    pot = params[0][1]
    vec = params[4][1]
    n=int(params[4][2:])
    if vec == 'g':
        d['particle'] = 'Cząstka o profilu gaussowskim wokół '+ str(n) + ' sigma=' + str(float(params[5]))
    elif vec == 'b':
        d['particle'] = 'Cząstka dodana do ' + str(n)
    if pot == 'l':
        d['potential'] = 'Potencjał: liniowy, V0='+str(V) + "*i"
    if pot == 'c':
        d['potential'] = 'Potencjał: stały, V0='+str(V)
    d["N"] = N
    d['dt'] = dt
    d['steps'] = s
    return d


def save_plot(filename, params):
    m = np.genfromtxt(filename, delimiter=",")
    m= np.flipud(m)
    fig, ax = plt.subplots()
    im = ax.imshow(m, cmap='viridis', origin='lower')
    plt.xticks(plt.xticks()[0][1::2], [int(x) for x in np.linspace(-params["N"]/2, params['N']/2, len(plt.xticks()[0][1::2]))])
    plt.yticks(plt.yticks()[0][1::2], np.linspace(0, params['steps']*params['dt'], len(plt.yticks()[0][1::2])))
    ax.set_xlabel("Miejsce w sieci")
    ax.set_ylabel("Czas   [j.u.]")
    
    #ax.set_title(params['potential']+"\n"+params['particle'])

    ax.tick_params(pad =5, bottom=False, left=False,)

    fig.colorbar(im, ax=ax, label='Prawdopodobieństwo', location="bottom", shrink=0.55, aspect = 30, pad = 0.18)
    
    fig.tight_layout()
    fig.savefig(filename[:-4]+".png", dpi=800)


def parse_files():
    for filename in os.listdir():
        if filename[:2] == "PM" and filename[-3:] == 'txt':
            params = params_dict(filename[:-4])    
            save_plot(filename, params)

parse_files()