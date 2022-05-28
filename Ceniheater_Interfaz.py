from tkinter import *
import numpy as np
from tkinter import ttk
import pyromat as pm
import math
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)

def mifuncion():
    def enviar_datos():
        #Lado jugo
        t_jugo_in = float(t_jugo_in_d.get())
        t_jugo_out_s = t_jugo_in+1
        pz = float(pz_d.get())
        bx = float(bx_d.get())

        mat_tubo = tubo_despl.get()
        if mat_tubo == 'Acero':
            c_t_tubo = 17
            rugosidad = 0.09
        elif mat_tubo == 'Cobre':
            c_t_tubo = 100
            rugosidad = 0.0015
        else:
            c_t_tubo = 64
            rugosidad = 0.24
        d_ext_s = d_ext_despl.get()
        if d_ext_s == 'Milimetros [mm]':
            d_ext = float(d_ext_d_entry.get())
        else:
             d_ext = float(d_ext_d_entry.get())*25.4

        d_int_s = d_int_despl.get()
        if d_int_s == 'Milimetros [mm]':
            d_int = float(d_int_d_entry.get())
        else:
             d_int = float(d_int_d_entry.get())*25.4

        long_tb_s = long_tb_despl.get()
        if long_tb_s == 'Metros [m]':
            long_tb = float(long_tb_d_entry.get())
        else:
             long_tb = float(long_tb_d_entry.get())*0.3048

        n_tubos = int(n_tubos_d.get())
        n_corazas = int(n_corazas_d.get())
        pasos_tubos = int(pasos_tubos_d.get())

        #Lado coraza
        t_cond_in = float(t_cond_d.get())
        t_cond_out_s = t_cond_in-1

        d_coraza_s = d_coraza_despl.get()
        if d_coraza_s == 'Diametro de coraza [m]':
            d_coraza = float(d_coraza_d_entry.get())
        else:
             d_coraza = float(d_coraza_d_entry.get())*25.4/1000

        area_crz = (3.141592654*d_coraza**2)/4
        area_efectiva = area_crz-(3.141592654*((d_ext/1000)**2)*n_tubos/4)
        long_caract = d_ext / 1000
        #long_caract = 2*math.sqrt(area_efectiva/(math.pi))

        #Incrustaciones
        h_operacion = float(h_operacion_d.get())
        r_ext_s = r_ext_despl.get()
        if r_ext_s == 'Agua por encima de 50°C':
            r_ext = 0.0002
        else:
            r_ext = 0.0001

        e=1
        tol=0.1

        while e > tol:
            t_mean_jugo = (t_jugo_in + t_jugo_out_s)/2
            #Calcula c_p del jugo
            cp_jugo = 4.1868-(bx*(0.0297-(0.000046*pz)))+(0.000075*bx*t_mean_jugo)
            #Calcula densidad del jugo
            d_jugo = 1005.3-(0.22556*t_mean_jugo)-(0.0024304*t_mean_jugo**2)+(3.7329*bx)+(0.01781937*bx**2)
            t_flujo = flujo_despl.get()
            if t_flujo == 'Flujo masico [t/h]':
                flujo = float(flujo_d.get())*1000/3600
            else:
                flujo = float(flujo_d.get())*d_jugo/3600
            #Conductividad termica jugo
            a = (5.466e-8*t_mean_jugo**2)-(1.176e-5*t_mean_jugo)-3.024e-3
            b = (-7.847e-6*t_mean_jugo**2)+(1.96e-3*t_mean_jugo)+0.563
            k_jugo = a*bx+b
            diff_term_j = k_jugo/(cp_jugo*d_jugo*1000)
            #Viscosidad
            aa = 0.85+(0.15*(pz/100))
            bb = bx*(0.962+0.038*(pz/100))
            cc = (30-t_mean_jugo)/(91+t_mean_jugo)
            dd = bb/(1900-(18*bb))
            v_j = (10**((22.46*dd)-0.14+(cc*(1.1+43.1*aa*dd**1.25))))*10**-3
            v_cin_j = v_j/d_jugo

            v_tubo = (((flujo/d_jugo)/(3.141592654*(d_int**2)/4)/(n_tubos/pasos_tubos)))*1000**2
            prandtl_jugo = v_j*cp_jugo*1000/k_jugo
            reynolds_jugo = v_tubo*(d_int/1000)/v_cin_j
            factor_friccion = (1/(-1.8*(np.log10((6.9/reynolds_jugo)+((rugosidad/d_int)/3.7)**1.11))))**2
            nusselt = ((factor_friccion/8)*(reynolds_jugo-1000)*prandtl_jugo)/(1+(12.7*((factor_friccion/8)**0.5)*(prandtl_jugo**(2/3)-1)))
            coef_conv_int = (k_jugo*nusselt)/(d_int/1000)

            #Coraza
            t_mean_cond = (t_cond_in + t_cond_out_s)/2
            #Incluir propiedades del agua
            H2O = pm.get('mp.H2O')
            #Densidad
            d_cond = H2O.d(T=(t_mean_cond+273),x=0)
            cp_cond = H2O.cp(T=(t_mean_cond+273),x=0)
            t_flujo_cond = flujo_cond_despl.get()
            if t_flujo == 'Flujo masico [t/h]':
                flujo_c = float(flujo_cond_d.get())*1000/3600
            else:
                flujo_c = float(flujo_cond_d.get())*d_cond/3600
            flujo_v_cond = flujo_c/d_cond

            k_agua = 0.5706+(0.001756*t_mean_cond)-(0.00000646*t_mean_cond**2)
            diff_term_cond = k_agua/(cp_cond*d_cond)
            if 0<t_mean_cond<=30:
                vsc_dinamica = 10**(-2.75-0.0141*t_mean_cond+(91.9e-6*t_mean_cond**2)-(311e-9*t_mean_cond**3))
            else:
                vsc_dinamica = (0.0168*d_cond*t_mean_cond**-0.88)/1000
            vsc_cinematica = vsc_dinamica/d_cond
            v_coraza = flujo_v_cond/area_efectiva
            prandtl_cond = vsc_dinamica*cp_cond*1000/k_agua
            reynolds_cond = v_coraza*long_caract/vsc_cinematica
            nusselt_cond = 0.024*(reynolds_cond**0.8)*(prandtl_cond**0.4)
            coef_conv_ext = (k_agua*nusselt_cond)/long_caract

            #Incrustaciones
            r_int = ((0.0035*h_operacion**0.8)*(1+(10.763/v_tubo**3)))/1000

            coef_global = 1/((1/coef_conv_int)+r_int+(np.log(d_ext/d_int)/(2*math.pi*c_t_tubo*long_tb))+r_ext+(1/coef_conv_ext))
            cc_jugo = flujo*cp_jugo
            cc_agua = flujo_c*cp_cond
            c_min = min([cc_jugo, cc_agua])
            c_max = max([cc_jugo, cc_agua])
            c_asterisco = c_min/c_max
            area_transfer_calor = (math.pi*(d_ext/1000)*long_tb*n_tubos)*n_corazas
            ntu = (coef_global*area_transfer_calor)/(c_min*1000)
            efectividad = 2*(1+c_asterisco+(math.sqrt(1+c_asterisco**2)*((1+math.exp(-ntu*math.sqrt(1+c_asterisco**2)))/(1-math.exp(-ntu*math.sqrt(1+c_asterisco**2))))))**-1
            q_max = c_min*(t_cond_in-t_jugo_in)
            q_real = efectividad*q_max
            t_out_jugo_calc = (q_real/cc_jugo)+t_jugo_in
            t_out_cond_calc = t_cond_in-(q_real/cc_agua)
            e=abs(t_out_cond_calc-t_cond_out_s)
            t_cond_out_s = t_out_cond_calc
            t_jugo_out_s = t_out_jugo_calc


        rpta_label = Label(ventana_STHX,text='Valores respuesta: ', bg='#ffffff', font=('Arial',12,'bold'))
        rpta_label.place(x=700, y=30)
        coef_global_label = Label(ventana_STHX,text='Coeficiente global de transferencia de calor: ' + str(round(coef_global[0],1)) + ' [W/m2°C]'+ str(round(c_min[0],1))+str(round(c_max[0],1)), bg='#ffffff', font=('Arial',12,'bold'))
        coef_global_label.place(x=700, y=60)
        area_transfer_calor_label = Label(ventana_STHX,text='Área de transferencia de calor: ' + str(round(area_transfer_calor,1)) + ' [m2]', bg='#ffffff', font=('Arial',12,'bold'))
        area_transfer_calor_label.place(x=700, y=90)
        ntu_label = Label(ventana_STHX,text='NTU: ' + str(round(ntu[0],3)), bg='#ffffff', font=('Arial',12,'bold'))
        ntu_label.place(x=700, y=120)
        v_tubo_label = Label(ventana_STHX,text='Velocidad del fluido en los tubos: ' + str(round(v_tubo[0],1)) + ' [m/s]', bg='#ffffff', font=('Arial',12,'bold'))
        v_tubo_label.place(x=700, y=150)
        efectividad_label = Label(ventana_STHX,text='Efectividad del intercambiador de calor: ' + str(round(efectividad[0]*100,2)) + ' [%]', bg='#ffffff', font=('Arial',12,'bold'))
        efectividad_label.place(x=700, y=180)
        q_max_label = Label(ventana_STHX,text='Transferencia maxima de calor posible: ' + str(round(q_max[0],1)) + ' [kW]', bg='#ffffff', font=('Arial',12,'bold'))
        q_max_label.place(x=700, y=210)
        q_real_label = Label(ventana_STHX,text='Transferencia de calor real: ' + str(round(q_real[0],1)) + ' [kW]', bg='#ffffff', font=('Arial',12,'bold'))
        q_real_label.place(x=700, y=240)
        t_in_jugo_label = Label(ventana_STHX,text='Temperatura de entrada del jugo: ' + str(round(t_jugo_in,1)) + ' [°C]', bg='#ffffff', font=('Arial',12,'bold'))
        t_in_jugo_label.place(x=700, y=270)
        t_out_jugo_calc_label = Label(ventana_STHX,text='Temperatura de salida del jugo: ' + str(round(t_out_jugo_calc[0],1)) + ' [°C]', bg='#ffffff', font=('Arial',12,'bold'))
        t_out_jugo_calc_label.place(x=700, y=300)
        t_in_cond_label = Label(ventana_STHX,text='Temperatura de entrada de los condensados: ' + str(round(t_cond_in,1)) + ' [°C]', bg='#ffffff', font=('Arial',12,'bold'))
        t_in_cond_label.place(x=700, y=330)
        t_out_cond_calc_label = Label(ventana_STHX,text='Temperatura de salida de los condensados: ' + str(round(t_out_cond_calc[0],1)) + ' [°C]', bg='#ffffff', font=('Arial',12,'bold'))
        t_out_cond_calc_label.place(x=700, y=360)

        fig = Figure(figsize = (5, 5), dpi = 100)
        x = np.array([0,100],dtype=object)
        y = np.array([t_out_jugo_calc, t_jugo_in],dtype=object)
        y2 = np.array([t_cond_in,t_out_cond_calc],dtype=object)
        plot1 = fig.add_subplot(111)
        plot1.plot(x,y,'<:',label='Jugo')
        plot1.plot(x,y2,'>:',label='Condensados')
        plot1.legend()
        plot1.grid(axis='y')
        #plt.title('Intercambiador de calor (Contraflujo)')
        #plt.xlabel("Porcentaje de superficie [%]")
        #plt.ylabel("Temperatura [°C]")
        #plt.show()
        canvas = FigureCanvasTkAgg(fig, master = ventana_STHX)
        canvas.draw()
        canvas.get_tk_widget().place(x=700,y=400)
        #toolbar = NavigationToolbar2Tk(canvas, ventana_STHX)
        #toolbar.update()
        #toolbar.get_tk_widget().place(x=700,y=800)
    ventana_STHX = Tk()
    ventana_STHX.iconbitmap('C:\\Users\\Dduqu\\Downloads\\Ceniheater_Python\\Ceniheater_Python\\cenicana.ico')
    ventana_STHX.geometry('1920x1080+0+0')
    ventana_STHX.title('Intercambiador de Coraza y Tubos')
    ventana_STHX.configure(bg='#ffffff')

    #Create a main frame
    main_frame = Frame(ventana_STHX)
    main_frame.pack(fill=BOTH, expand=1)
    main_frame.configure(bg='#ffffff')
    #Create a Canvas
    my_canvas = Canvas(main_frame)
    my_canvas.pack(side=LEFT, fill=BOTH, expand=1)
    #Add a scrollbar to the canvas
    my_scrollbar = ttk.Scrollbar(main_frame, orient=VERTICAL, command=my_canvas.yview)
    my_scrollbar.pack(side=RIGHT, fill=Y)
    #Configure the Canvas
    my_canvas.configure(bg='#ffffff',yscrollcommand=my_scrollbar.set)
    my_canvas.bind('<Configure>', lambda e: my_canvas.configure(scrollregion = my_canvas.bbox('all')))
    #Create another frame inside the Canvas
    second_frame = Frame(my_canvas)
    second_frame.configure(bg='#ffffff')
    #Add that new frame to a window in the canvas
    my_canvas.create_window((0,0), window=second_frame, anchor='nw')


    lbl = Label(ventana_STHX,text = 'Bienvenido, al programa de calculo de intercambiadores de Calor Coraza-Tubos Condensados/Jugo', font=('Arial',12))
    lbl.configure(bg='#ffffff')
    lbl.pack()
    #Inicialización de variables
    t_jugo_in = 0
    t_jugo_out = 0
    t_mean_jugo = 0
    pz = 0
    bx = 0
    cp_jugo = 0
    d = 0
    flujo = 0
    d_ext = 0
    d_int = 0
    long_tb = 0
    n_tubos = 0
    n_corazas = 0
    pasos_tubos = 0
    t_cond = 0
    flujo_cond = 0
    d_coraza = 0
    h_operacion = 0
    #Ingresa datos
    #Datos jugos
    info_label = Label(ventana_STHX,text='Datos Jugo - Interior tubos:', bg='#ffffff', font=('Arial',12,'bold'))
    info_label.place(x=10, y=30)
    t_jugo_in_label = Label(ventana_STHX,text='Temperatura de entrada del jugo [°C] :', bg='#ffffff', font=('Arial',12))
    t_jugo_in_label.place(x=10, y=60)
    pz_label = Label(ventana_STHX,text='Pureza del jugo [%] :', bg='#ffffff', font=('Arial',12))
    pz_label.place(x=10, y=90)
    bx_label = Label(ventana_STHX,text='Brix del jugo [°] :', bg='#ffffff', font=('Arial',12))
    bx_label.place(x=10, y=120)
    flujo_label = Label(ventana_STHX,text='Seleccione el tipo de flujo e ingreselo:', bg='#ffffff', font=('Arial',12))
    flujo_label.place(x=10, y=150)
    flujo_despl = ttk.Combobox(ventana_STHX)
    flujo_despl.place(x=10,y=180)
    flujo_despl['values'] = ('Flujo volumetrico [m3/h]', 'Flujo masico [t/h]')
    flujo_despl.current(0)
    #Datos de los tubos
    tubo_label = Label(ventana_STHX,text='Seleccione el material de los tubos:', bg='#ffffff', font=('Arial',12))
    tubo_label.place(x=10, y=210)
    tubo_despl = ttk.Combobox(ventana_STHX)
    tubo_despl.place(x=270,y=212)
    tubo_despl['values'] = ('Acero', 'Cobre', 'Niquel')
    tubo_despl.current(0)
    d_ext_label = Label(ventana_STHX,text='Seleccione la unidad del diametro externo de tubos e ingrese el valor:', bg='#ffffff', font=('Arial',12))
    d_ext_label.place(x=10, y=240)
    d_ext_despl = ttk.Combobox(ventana_STHX)
    d_ext_despl.place(x=10,y=270)
    d_ext_despl['values'] = ('Milimetros [mm]', 'Pulgada [in]')
    d_ext_despl.current(0)
    d_int_label = Label(ventana_STHX,text='Seleccione la unidad del diametro interno de tubos e ingrese el valor:', bg='#ffffff', font=('Arial',12))
    d_int_label.place(x=10, y=300)
    d_int_despl = ttk.Combobox(ventana_STHX)
    d_int_despl.place(x=10,y=330)
    d_int_despl['values'] = ('Milimetros [mm]', 'Pulgada [in]')
    d_int_despl.current(0)
    long_tb_label = Label(ventana_STHX,text='Seleccione la unidad de longitud de los tubos e ingrese el valor:', bg='#ffffff', font=('Arial',12))
    long_tb_label.place(x=10, y=360)
    long_tb_despl = ttk.Combobox(ventana_STHX)
    long_tb_despl.place(x=10,y=390)
    long_tb_despl['values'] = ('Metros [m]', 'Píes [ft]')
    long_tb_despl.current(0)
    n_tubos_label = Label(ventana_STHX,text='Cantidad total de tubos:', bg='#ffffff', font=('Arial',12))
    n_tubos_label.place(x=10, y=420)
    n_corazas_label = Label(ventana_STHX,text='Cantidad de corazas en serie:', bg='#ffffff', font=('Arial',12))
    n_corazas_label.place(x=10, y=450)
    pasos_tubos_label = Label(ventana_STHX,text='Pasos de tubos por coraza:', bg='#ffffff', font=('Arial',12))
    pasos_tubos_label.place(x=10, y=480)
    #Datos Condensados
    info1_label = Label(ventana_STHX,text='Datos Condensados - Lado coraza:', bg='#ffffff', font=('Arial',12,'bold'))
    info1_label.place(x=10, y=510)
    t_cond_in_label = Label(ventana_STHX,text='Temperatura de entrada de los condensados [°C]:', bg='#ffffff', font=('Arial',12))
    t_cond_in_label.place(x=10, y=540)
    flujo_cond_label = Label(ventana_STHX,text='Seleccione el tipo de flujo e ingreselo:', bg='#ffffff', font=('Arial',12))
    flujo_cond_label.place(x=10, y=570)
    flujo_cond_despl = ttk.Combobox(ventana_STHX)
    flujo_cond_despl.place(x=10,y=600)
    flujo_cond_despl['values'] = ('Flujo volumetrico [m3/h]', 'Flujo masico [t/h]')
    flujo_cond_despl.current(0)
    d_coraza_label = Label(ventana_STHX,text='Diametro de la coraza:', bg='#ffffff', font=('Arial',12))
    d_coraza_label.place(x=10, y=630)
    d_coraza_despl = ttk.Combobox(ventana_STHX)
    d_coraza_despl.place(x=10,y=660)
    d_coraza_despl['values'] = ('Diametro de coraza [m]', 'Diametro de coraza [in]')
    d_coraza_despl.current(0)
    #Factores de incrustación
    info2_label = Label(ventana_STHX,text='Datos Incrustaciones:', bg='#ffffff', font=('Arial',12,'bold'))
    info2_label.place(x=10, y=690)
    h_operacion_label = Label(ventana_STHX,text='Horas de operación:', bg='#ffffff', font=('Arial',12))
    h_operacion_label.place(x=10, y=720)
    r_ext_label = Label(ventana_STHX,text='Seleccione el fluido externo (Lado coraza):', bg='#ffffff', font=('Arial',12))
    r_ext_label.place(x=10, y=750)
    r_ext_despl = ttk.Combobox(ventana_STHX)
    r_ext_despl.place(x=320,y=752)
    r_ext_despl['values'] = ('Agua por debajo de 50°C', 'Agua por encima de 50°C', 'Vapor de Agua')
    r_ext_despl.current(1)

    #Capturar datos de entrada
    #Temperatura del jugo
    t_jugo_in_d = DoubleVar(ventana_STHX)
    t_jugo_in_d_entry = Entry(ventana_STHX,textvariable=t_jugo_in_d, width='15')
    t_jugo_in_d_entry.place(x=280,y=65)
    #Pureza
    pz_d = DoubleVar(ventana_STHX)
    pz_d_entry = Entry(ventana_STHX,textvariable=pz_d, width='15')
    pz_d_entry.place(x=160,y=95)
    #Brix
    bx_d = DoubleVar(ventana_STHX)
    bx_d_entry = Entry(ventana_STHX,textvariable=bx_d, width='15')
    bx_d_entry.place(x=140,y=125)
    #Flujo
    flujo_d = DoubleVar(ventana_STHX)
    flujo_d_entry = Entry(ventana_STHX,textvariable=flujo_d, width='15')
    flujo_d_entry.place(x=160,y=180)
    #Diametro Externo
    d_ext_d = DoubleVar(ventana_STHX)
    d_ext_d_entry = Entry(ventana_STHX,textvariable=d_ext_d, width='15')
    d_ext_d_entry.place(x=160,y=270)
    #Diametro Interno
    d_int_d = DoubleVar(ventana_STHX)
    d_int_d_entry = Entry(ventana_STHX,textvariable=d_int_d, width='15')
    d_int_d_entry.place(x=160,y=330)
    #Longitud de Tubos
    long_tb_d = DoubleVar(ventana_STHX)
    long_tb_d_entry = Entry(ventana_STHX,textvariable=long_tb_d, width='15')
    long_tb_d_entry.place(x=160,y=390)
    #Cantidad de tubos
    n_tubos_d = DoubleVar(ventana_STHX)
    n_tubos_d_entry = Entry(ventana_STHX,textvariable=n_tubos_d, width='15')
    n_tubos_d_entry.place(x=185,y=420)
    #Cantidad de corazas en serie
    n_corazas_d = DoubleVar(ventana_STHX)
    n_corazas_d_entry = Entry(ventana_STHX,textvariable=n_corazas_d, width='15')
    n_corazas_d_entry.place(x=230,y=450)
    #Pasos de tubos por Coraza
    pasos_tubos_d = DoubleVar(ventana_STHX)
    pasos_tubos_d_entry = Entry(ventana_STHX,textvariable=pasos_tubos_d, width='15')
    pasos_tubos_d_entry.place(x=215,y=480)
    #Temperatura entrada condensados
    t_cond_d = DoubleVar(ventana_STHX)
    t_cond_d_entry = Entry(ventana_STHX,textvariable=t_cond_d, width='15')
    t_cond_d_entry.place(x=365,y=542)
    #Flujo de condensados
    flujo_cond_d = DoubleVar(ventana_STHX)
    flujo_cond_d_entry = Entry(ventana_STHX,textvariable=flujo_cond_d, width='15')
    flujo_cond_d_entry.place(x=160,y=600)
    #Diametro de la coraza
    d_coraza_d = DoubleVar(ventana_STHX)
    d_coraza_d_entry = Entry(ventana_STHX,textvariable=d_coraza_d, width='15')
    d_coraza_d_entry.place(x=160,y=660)
    #Horas de operación
    h_operacion_d = DoubleVar(ventana_STHX)
    h_operacion_d_entry = Entry(ventana_STHX,textvariable=h_operacion_d, width='15')
    h_operacion_d_entry.place(x=160,y=720)
    #Botón Calcular
    sbmit_bttn = Button(ventana_STHX, text='Calcular', command = enviar_datos, width='30', height='2', bg='#ffa600',font=('Arial',12))
    sbmit_bttn.place(x=25,y=800)

    ventana_STHX.mainloop()

def mifuncion2():
    ventana_PHX = Tk()
    ventana_PHX.iconbitmap('C:\\Users\\Dduqu\\Downloads\\Ceniheater_Python\\Ceniheater_Python\\cenicana.ico')
    ventana_PHX.geometry('900x900')
    ventana_PHX.title('Intercambiador de Placas')
    ventana_PHX.configure(bg='#ffffff')

    #Create a main frame
    #main_frame = Frame(ventana_PHX)
    #main_frame.pack(fill=BOTH, expand=1)
    #main_frame.configure(bg='#ffffff')
    #Create a Canvas
    #my_canvas = Canvas(main_frame)
    #my_canvas.pack(side=LEFT, fill=BOTH, expand=1)
    #Add a scrollbar to the canvas
    #my_scrollbar = ttk.Scrollbar(main_frame, orient=VERTICAL, command=my_canvas.yview)
    #my_scrollbar.pack(side=RIGHT, fill=Y)
    #Configure the Canvas
    #my_canvas.configure(bg='#ffffff',yscrollcommand=my_scrollbar.set)
    #my_canvas.bind('<Configure>', lambda e: my_canvas.configure(scrollregion = my_canvas.bbox('all')))
    #Create another frame inside the Canvas
    #second_frame = Frame(my_canvas)
    #second_frame.configure(bg='#ffffff')
    #Add that new frame to a window in the canvas
    #my_canvas.create_window((0,0), window=second_frame, anchor='nw')

    #for i in range(100):
    #    Button(second_frame, text=f'Button {i} Yo!').grid(row=i, column=0, pady=10, padx=10)
def mifuncion3():
    ventana_AHX = Tk()
    ventana_AHX.iconbitmap('C:\\Users\\Dduqu\\Downloads\\Ceniheater_Python\\Ceniheater_Python\\cenicana.ico')
    ventana_AHX.geometry('900x900')
    ventana_AHX.title('Intercambiador con Aletas')
    ventana_AHX.configure(bg='#ffffff')
def mifuncion4():
    def enviar_datos2():
        #Lado jugo
        t_jugo_in = float(t_jugo_in_d.get())
        t_jugo_out = float(t_jugo_out_d.get())
        pz = float(pz_d.get())
        bx = float(bx_d.get())

        #mat_tubo = tubo_despl.get()
        #if mat_tubo == 'Acero':
        #    c_t_tubo = 17
        #    rugosidad = 0.09
        #elif mat_tubo == 'Cobre':
        #    c_t_tubo = 100
        #    rugosidad = 0.0015
        #else:
        #    c_t_tubo = 64
        #    rugosidad = 0.24
        d_ext_s = d_ext_despl.get()
        if d_ext_s == 'Milimetros [mm]':
            d_ext = float(d_ext_d_entry.get())
        else:
             d_ext = float(d_ext_d_entry.get())*25.4

        d_int_s = d_int_despl.get()
        if d_int_s == 'Milimetros [mm]':
            d_int = float(d_int_d_entry.get())
        else:
             d_int = float(d_int_d_entry.get())*25.4

        long_tb_s = long_tb_despl.get()
        if long_tb_s == 'Metros [m]':
            long_tb = float(long_tb_d_entry.get())
        else:
             long_tb = float(long_tb_d_entry.get())*0.3048

        n_tubos = int(n_tubos_d.get())
        n_corazas = int(n_corazas_d.get())
        pasos_tubos = int(pasos_tubos_d.get())

        #Lado coraza
        t_cond_in = float(t_cond_d.get())
        t_cond_out = float(t_cond_out_d.get())

        #d_coraza_s = d_coraza_despl.get()
        #if d_coraza_s == 'Diametro de coraza [m]':
        #    d_coraza = float(d_coraza_d_entry.get())
        #else:
        #     d_coraza = float(d_coraza_d_entry.get())*25.4/1000

        #area_crz = (3.141592654*d_coraza**2)/4
        #area_efectiva = area_crz-(3.141592654*((d_ext/1000)**2)*n_tubos/4)
        #long_caract = d_ext / 1000
        #long_caract = 2*math.sqrt(area_efectiva/(math.pi))

        #Incrustaciones
        #h_operacion = float(h_operacion_d.get())
        #r_ext_s = r_ext_despl.get()
        #if r_ext_s == 'Agua por encima de 50°C':
        #    r_ext = 0.0002
        #else:
        #    r_ext = 0.0001

        t_mean_jugo = (t_jugo_in + t_jugo_out)/2
        #Calcula c_p del jugo
        cp_jugo = 4.1868-(bx*(0.0297-(0.000046*pz)))+(0.000075*bx*t_mean_jugo)
        #Calcula densidad del jugo
        d_jugo = 1005.3-(0.22556*t_mean_jugo)-(0.0024304*t_mean_jugo**2)+(3.7329*bx)+(0.01781937*bx**2)
        t_flujo = flujo_despl.get()
        if t_flujo == 'Flujo masico [t/h]':
            flujo = float(flujo_d.get())*1000/3600
        else:
            flujo = float(flujo_d.get())*d_jugo/3600
        #Conductividad termica jugo
        #a = (5.466e-8*t_mean_jugo**2)-(1.176e-5*t_mean_jugo)-3.024e-3
        #b = (-7.847e-6*t_mean_jugo**2)+(1.96e-3*t_mean_jugo)+0.563
        #k_jugo = a*bx+b
        #diff_term_j = k_jugo/(cp_jugo*d_jugo*1000)
        #Viscosidad
        #aa = 0.85+(0.15*(pz/100))
        #bb = bx*(0.962+0.038*(pz/100))
        #cc = (30-t_mean_jugo)/(91+t_mean_jugo)
        #dd = bb/(1900-(18*bb))
        #v_j = (10**((22.46*dd)-0.14+(cc*(1.1+43.1*aa*dd**1.25))))*10**-3
        #v_cin_j = v_j/d_jugo

        v_tubo = (((flujo/d_jugo)/(3.141592654*(d_int**2)/4)/(n_tubos/pasos_tubos)))*1000**2
        #prandtl_jugo = v_j*cp_jugo*1000/k_jugo
        #reynolds_jugo = v_tubo*(d_int/1000)/v_cin_j
        #factor_friccion = (1/(-1.8*(np.log10((6.9/reynolds_jugo)+((rugosidad/d_int)/3.7)**1.11))))**2
        #nusselt = ((factor_friccion/8)*(reynolds_jugo-1000)*prandtl_jugo)/(1+(12.7*((factor_friccion/8)**0.5)*(prandtl_jugo**(2/3)-1)))
        #coef_conv_int = (k_jugo*nusselt)/(d_int/1000)

        #Coraza
        t_mean_cond = (t_cond_in + t_cond_out)/2
        #Incluir propiedades del agua
        H2O = pm.get('mp.H2O')
        #Densidad
        d_cond = H2O.d(T=(t_mean_cond+273),x=0)
        cp_cond = H2O.cp(T=(t_mean_cond+273),x=0)
        #t_flujo_cond = flujo_cond_despl.get()
        #if t_flujo == 'Flujo masico [t/h]':
        #        flujo_c = float(flujo_cond_d.get())*1000/3600
        #    else:
        #        flujo_c = float(flujo_cond_d.get())*d_cond/3600
        #    flujo_v_cond = flujo_c/d_cond

        #    k_agua = 0.5706+(0.001756*t_mean_cond)-(0.00000646*t_mean_cond**2)
        #    diff_term_cond = k_agua/(cp_cond*d_cond)
        #    if 0<t_mean_cond<=30:
        #        vsc_dinamica = 10**(-2.75-0.0141*t_mean_cond+(91.9e-6*t_mean_cond**2)-(311e-9*t_mean_cond**3))
        #    else:
        #        vsc_dinamica = (0.0168*d_cond*t_mean_cond**-0.88)/1000
        #    vsc_cinematica = vsc_dinamica/d_cond
        #    v_coraza = flujo_v_cond/area_efectiva
        #    prandtl_cond = vsc_dinamica*cp_cond*1000/k_agua
        #    reynolds_cond = v_coraza*long_caract/vsc_cinematica
        #    nusselt_cond = 0.024*(reynolds_cond**0.8)*(prandtl_cond**0.4)
        #    coef_conv_ext = (k_agua*nusselt_cond)/long_caract

            #Incrustaciones
        #   r_int = ((0.0035*h_operacion**0.8)*(1+(10.763/v_tubo**3)))/1000

        #    coef_global = 1/((1/coef_conv_int)+r_int+(np.log(d_ext/d_int)/(2*math.pi*c_t_tubo*long_tb))+r_ext+(1/coef_conv_ext))
        cc_jugo = flujo*cp_jugo
        q_real = cc_jugo*(t_jugo_out-t_jugo_in)
        flujo_c = (q_real/(H2O.h(T=t_cond_in+273.15, x=0)-H2O.h(T=t_cond_out+273.15, x=0)))*3.6
        cc_agua = flujo_c*cp_cond
        c_min = min([cc_jugo, cc_agua])
        c_max = max([cc_jugo, cc_agua])
        c_asterisco = c_min/c_max
        area_transfer_calor = (math.pi*(d_ext/1000)*long_tb*n_tubos)*n_corazas
        q_max = c_min*(t_cond_in-t_jugo_in)
        efectividad = q_real/q_max
        dta=(t_cond_in-t_jugo_in)
        dtb=(t_cond_out-t_jugo_out)
        dtln = (dta-dtb)/(np.log(dta/dtb))
        coef_global = (q_real*1000) /(area_transfer_calor*dtln)
        ntu = (coef_global*area_transfer_calor)/(c_min*1000)


        rpta_label = Label(ventana_STHX2,text='Valores respuesta: ', bg='#ffffff', font=('Arial',12,'bold'))
        rpta_label.place(x=700, y=30)
        coef_global_label = Label(ventana_STHX2,text='Coeficiente global de transferencia de calor: ' + str(round(coef_global,2)) + ' [W/m2°C]', bg='#ffffff', font=('Arial',12,'bold'))
        coef_global_label.place(x=700, y=60)
        area_transfer_calor_label = Label(ventana_STHX2,text='Área de transferencia de calor: ' + str(round(area_transfer_calor,1)) + ' [m2]', bg='#ffffff', font=('Arial',12,'bold'))
        area_transfer_calor_label.place(x=700, y=90)
        ntu_label = Label(ventana_STHX2,text='Requerimiento de condensados: ' + str(round(flujo_c[0],1)) + '[t/h]', bg='#ffffff', font=('Arial',12,'bold'))
        ntu_label.place(x=700, y=120)
        v_tubo_label = Label(ventana_STHX2,text='Velocidad del fluido en los tubos: ' + str(round(v_tubo,1)) + ' [m/s]', bg='#ffffff', font=('Arial',12,'bold'))
        v_tubo_label.place(x=700, y=150)
        efectividad_label = Label(ventana_STHX2,text='Efectividad del intercambiador de calor: ' + str(round(efectividad*100,2)) + ' [%]', bg='#ffffff', font=('Arial',12,'bold'))
        efectividad_label.place(x=700, y=180)
        q_max_label = Label(ventana_STHX2,text='Transferencia maxima de calor posible: ' + str(round(q_max,1)) + ' [kW]', bg='#ffffff', font=('Arial',12,'bold'))
        q_max_label.place(x=700, y=210)
        q_real_label = Label(ventana_STHX2,text='Transferencia de calor real: ' + str(round(q_real,1)) + ' [kW]', bg='#ffffff', font=('Arial',12,'bold'))
        q_real_label.place(x=700, y=240)
        c_label = Label(ventana_STHX2,text='C_min & C_max: ' + str(c_min[0]) + str(c_max[0]), bg='#ffffff', font=('Arial',12,'bold'))
        c_label.place(x=700, y=270)
        #t_in_jugo_label = Label(ventana_STHX,text='Temperatura de entrada del jugo: ' + str(round(t_jugo_in,1)) + ' [°C]', bg='#ffffff', font=('Arial',12,'bold'))
        #t_in_jugo_label.place(x=700, y=270)
        #t_out_jugo_calc_label = Label(ventana_STHX,text='Temperatura de salida del jugo: ' + str(round(t_out_jugo_calc[0],1)) + ' [°C]', bg='#ffffff', font=('Arial',12,'bold'))
        #t_out_jugo_calc_label.place(x=700, y=300)
        #t_in_cond_label = Label(ventana_STHX,text='Temperatura de entrada de los condensados: ' + str(round(t_cond_in,1)) + ' [°C]', bg='#ffffff', font=('Arial',12,'bold'))
        #t_in_cond_label.place(x=700, y=330)
        #t_out_cond_calc_label = Label(ventana_STHX,text='Temperatura de salida de los condensados: ' + str(round(t_out_cond_calc[0],1)) + ' [°C]', bg='#ffffff', font=('Arial',12,'bold'))
        #t_out_cond_calc_label.place(x=700, y=360)

        fig = Figure(figsize = (5, 5), dpi = 100)
        x = np.array([0,100],dtype=object)
        y = np.array([t_jugo_out, t_jugo_in],dtype=object)
        y2 = np.array([t_cond_in,t_cond_out],dtype=object)
        plot1 = fig.add_subplot(111)
        plot1.plot(x,y,'<:',label='Jugo')
        plot1.plot(x,y2,'>:',label='Condensados')
        plot1.legend()
        plot1.grid(axis='y')
        #plt.title('Intercambiador de calor (Contraflujo)')
        #plt.xlabel("Porcentaje de superficie [%]")
        #plt.ylabel("Temperatura [°C]")
        #plt.show()
        canvas = FigureCanvasTkAgg(fig, master = ventana_STHX2)
        canvas.draw()
        canvas.get_tk_widget().place(x=700,y=400)
        #toolbar = NavigationToolbar2Tk(canvas, ventana_STHX)
        #toolbar.update()
        #toolbar.get_tk_widget().place(x=700,y=800)
    ventana_STHX2 = Tk()
    ventana_STHX2.iconbitmap('C:\\Users\\Dduqu\\Downloads\\Ceniheater_Python\\Ceniheater_Python\\cenicana.ico')
    ventana_STHX2.geometry('1920x1080+0+0')
    ventana_STHX2.title('Intercambiador de Coraza y Tubos: Analisis de Temperatura a flujo')
    ventana_STHX2.configure(bg='#ffffff')

    lbl = Label(ventana_STHX2,text = 'Bienvenido, al programa de calculo de intercambiadores de Calor Coraza-Tubos Condensados/Jugo', font=('Arial',12))
    lbl.configure(bg='#ffffff')
    lbl.pack()
    #Inicialización de variables
    t_jugo_in = 0
    t_jugo_out = 0
    t_mean_jugo = 0
    pz = 0
    bx = 0
    cp_jugo = 0
    d = 0
    flujo = 0
    d_ext = 0
    d_int = 0
    long_tb = 0
    n_tubos = 0
    n_corazas = 0
    pasos_tubos = 0
    t_cond = 0
    flujo_cond = 0
    d_coraza = 0
    h_operacion = 0
    #Ingresa datos
    #Datos jugos
    info_label = Label(ventana_STHX2,text='Datos Jugo - Interior tubos:', bg='#ffffff', font=('Arial',12,'bold'))
    info_label.place(x=10, y=30)
    t_jugo_in_label = Label(ventana_STHX2,text='Temperatura de entrada del jugo [°C] :', bg='#ffffff', font=('Arial',12))
    t_jugo_in_label.place(x=10, y=60)
    pz_label = Label(ventana_STHX2,text='Pureza del jugo [%] :', bg='#ffffff', font=('Arial',12))
    pz_label.place(x=10, y=90)
    bx_label = Label(ventana_STHX2,text='Brix del jugo [°] :', bg='#ffffff', font=('Arial',12))
    bx_label.place(x=10, y=120)
    flujo_label = Label(ventana_STHX2,text='Seleccione el tipo de flujo e ingreselo:', bg='#ffffff', font=('Arial',12))
    flujo_label.place(x=10, y=150)
    flujo_despl = ttk.Combobox(ventana_STHX2)
    flujo_despl.place(x=10,y=180)
    flujo_despl['values'] = ('Flujo volumetrico [m3/h]', 'Flujo masico [t/h]')
    flujo_despl.current(0)
    #Datos de los tubos
    tubo_label = Label(ventana_STHX2,text='Seleccione el material de los tubos:', bg='#ffffff', font=('Arial',12))
    tubo_label.place(x=10, y=210)
    tubo_despl = ttk.Combobox(ventana_STHX2)
    tubo_despl.place(x=270,y=212)
    tubo_despl['values'] = ('Acero', 'Cobre', 'Niquel')
    tubo_despl.current(0)
    d_ext_label = Label(ventana_STHX2,text='Seleccione la unidad del diametro externo de tubos e ingrese el valor:', bg='#ffffff', font=('Arial',12))
    d_ext_label.place(x=10, y=240)
    d_ext_despl = ttk.Combobox(ventana_STHX2)
    d_ext_despl.place(x=10,y=270)
    d_ext_despl['values'] = ('Milimetros [mm]', 'Pulgada [in]')
    d_ext_despl.current(0)
    d_int_label = Label(ventana_STHX2,text='Seleccione la unidad del diametro interno de tubos e ingrese el valor:', bg='#ffffff', font=('Arial',12))
    d_int_label.place(x=10, y=300)
    d_int_despl = ttk.Combobox(ventana_STHX2)
    d_int_despl.place(x=10,y=330)
    d_int_despl['values'] = ('Milimetros [mm]', 'Pulgada [in]')
    d_int_despl.current(0)
    long_tb_label = Label(ventana_STHX2,text='Seleccione la unidad de longitud de los tubos e ingrese el valor:', bg='#ffffff', font=('Arial',12))
    long_tb_label.place(x=10, y=360)
    long_tb_despl = ttk.Combobox(ventana_STHX2)
    long_tb_despl.place(x=10,y=390)
    long_tb_despl['values'] = ('Metros [m]', 'Píes [ft]')
    long_tb_despl.current(0)
    n_tubos_label = Label(ventana_STHX2,text='Cantidad total de tubos:', bg='#ffffff', font=('Arial',12))
    n_tubos_label.place(x=10, y=420)
    n_corazas_label = Label(ventana_STHX2,text='Cantidad de corazas en serie:', bg='#ffffff', font=('Arial',12))
    n_corazas_label.place(x=10, y=450)
    pasos_tubos_label = Label(ventana_STHX2,text='Pasos de tubos por coraza:', bg='#ffffff', font=('Arial',12))
    pasos_tubos_label.place(x=10, y=480)
    t_jugo_out_label = Label(ventana_STHX2,text='Temperatura de salida del jugo [°C] :', bg='#ffffff', font=('Arial',12))
    t_jugo_out_label.place(x=10, y=510)
    #Datos Condensados
    info1_label = Label(ventana_STHX2,text='Datos Condensados - Lado coraza:', bg='#ffffff', font=('Arial',12,'bold'))
    info1_label.place(x=10, y=540)
    t_cond_in_label = Label(ventana_STHX2,text='Temperatura de entrada de los condensados [°C]:', bg='#ffffff', font=('Arial',12))
    t_cond_in_label.place(x=10, y=570)
    t_cond_out_label = Label(ventana_STHX2,text='Temperatura de salida de los condensados [°C]:', bg='#ffffff', font=('Arial',12))
    t_cond_out_label.place(x=10, y=600)
    #flujo_cond_label = Label(ventana_STHX2,text='Seleccione el tipo de flujo e ingreselo:', bg='#ffffff', font=('Arial',12))
    #flujo_cond_label.place(x=10, y=570)
    #flujo_cond_despl = ttk.Combobox(ventana_STHX2)
    #flujo_cond_despl.place(x=10,y=600)
    #flujo_cond_despl['values'] = ('Flujo volumetrico [m3/h]', 'Flujo masico [t/h]')
    #flujo_cond_despl.current(0)
    d_coraza_label = Label(ventana_STHX2,text='Diametro de la coraza:', bg='#ffffff', font=('Arial',12))
    d_coraza_label.place(x=10, y=630)
    d_coraza_despl = ttk.Combobox(ventana_STHX2)
    d_coraza_despl.place(x=10,y=660)
    d_coraza_despl['values'] = ('Diametro de coraza [m]', 'Diametro de coraza [in]')
    d_coraza_despl.current(0)
    #Factores de incrustación
    info2_label = Label(ventana_STHX2,text='Datos Incrustaciones:', bg='#ffffff', font=('Arial',12,'bold'))
    info2_label.place(x=10, y=690)
    h_operacion_label = Label(ventana_STHX2,text='Horas de operación:', bg='#ffffff', font=('Arial',12))
    h_operacion_label.place(x=10, y=720)
    r_ext_label = Label(ventana_STHX2,text='Seleccione el fluido externo (Lado coraza):', bg='#ffffff', font=('Arial',12))
    r_ext_label.place(x=10, y=750)
    r_ext_despl = ttk.Combobox(ventana_STHX2)
    r_ext_despl.place(x=320,y=752)
    r_ext_despl['values'] = ('Agua por debajo de 50°C', 'Agua por encima de 50°C', 'Vapor de Agua')
    r_ext_despl.current(1)

    #Capturar datos de entrada
    #Temperatura del jugo
    t_jugo_in_d = DoubleVar(ventana_STHX2)
    t_jugo_in_d_entry = Entry(ventana_STHX2,textvariable=t_jugo_in_d, width='15')
    t_jugo_in_d_entry.place(x=280,y=65)
    t_jugo_out_d = DoubleVar(ventana_STHX2)
    t_jugo_out_d_entry = Entry(ventana_STHX2,textvariable=t_jugo_out_d, width='15')
    t_jugo_out_d_entry.place(x=280,y=510)
    #Pureza
    pz_d = DoubleVar(ventana_STHX2)
    pz_d_entry = Entry(ventana_STHX2,textvariable=pz_d, width='15')
    pz_d_entry.place(x=160,y=95)
    #Brix
    bx_d = DoubleVar(ventana_STHX2)
    bx_d_entry = Entry(ventana_STHX2,textvariable=bx_d, width='15')
    bx_d_entry.place(x=140,y=125)
    #Flujo
    flujo_d = DoubleVar(ventana_STHX2)
    flujo_d_entry = Entry(ventana_STHX2,textvariable=flujo_d, width='15')
    flujo_d_entry.place(x=160,y=180)
    #Diametro Externo
    d_ext_d = DoubleVar(ventana_STHX2)
    d_ext_d_entry = Entry(ventana_STHX2,textvariable=d_ext_d, width='15')
    d_ext_d_entry.place(x=160,y=270)
    #Diametro Interno
    d_int_d = DoubleVar(ventana_STHX2)
    d_int_d_entry = Entry(ventana_STHX2,textvariable=d_int_d, width='15')
    d_int_d_entry.place(x=160,y=330)
    #Longitud de Tubos
    long_tb_d = DoubleVar(ventana_STHX2)
    long_tb_d_entry = Entry(ventana_STHX2,textvariable=long_tb_d, width='15')
    long_tb_d_entry.place(x=160,y=390)
    #Cantidad de tubos
    n_tubos_d = DoubleVar(ventana_STHX2)
    n_tubos_d_entry = Entry(ventana_STHX2,textvariable=n_tubos_d, width='15')
    n_tubos_d_entry.place(x=185,y=420)
    #Cantidad de corazas en serie
    n_corazas_d = DoubleVar(ventana_STHX2)
    n_corazas_d_entry = Entry(ventana_STHX2,textvariable=n_corazas_d, width='15')
    n_corazas_d_entry.place(x=230,y=450)
    #Pasos de tubos por Coraza
    pasos_tubos_d = DoubleVar(ventana_STHX2)
    pasos_tubos_d_entry = Entry(ventana_STHX2,textvariable=pasos_tubos_d, width='15')
    pasos_tubos_d_entry.place(x=215,y=480)
    #Temperatura entrada condensados
    t_cond_d = DoubleVar(ventana_STHX2)
    t_cond_d_entry = Entry(ventana_STHX2,textvariable=t_cond_d, width='15')
    t_cond_d_entry.place(x=365,y=572)
    #Temperatura salida condensados
    t_cond_out_d = DoubleVar(ventana_STHX2)
    t_cond_out_d_entry = Entry(ventana_STHX2,textvariable=t_cond_out_d, width='15')
    t_cond_out_d_entry.place(x=365,y=602)
    #Diametro de la coraza
    d_coraza_d = DoubleVar(ventana_STHX2)
    d_coraza_d_entry = Entry(ventana_STHX2,textvariable=d_coraza_d, width='15')
    d_coraza_d_entry.place(x=160,y=660)
    #Horas de operación
    h_operacion_d = DoubleVar(ventana_STHX2)
    h_operacion_d_entry = Entry(ventana_STHX2,textvariable=h_operacion_d, width='15')
    h_operacion_d_entry.place(x=160,y=720)
    #Botón Calcular
    sbmit_bttn = Button(ventana_STHX2, text='Calcular', command = enviar_datos2, width='30', height='2', bg='#ffa600',font=('Arial',12))
    sbmit_bttn.place(x=25,y=800)

    ventana_STHX2.mainloop()

ventana = Tk()
ventana.iconbitmap('C:\\Users\\Dduqu\\Downloads\\Ceniheater_Python\\Ceniheater_Python\\cenicana.ico')
ventana.geometry('900x900+0+0')
#ventana.resizable(False,False)
ventana.title('CeniHetaer')
ventana.configure(bg='#ffffff')

lbl = Label(ventana,text = 'Bienvenido al programa Ceniheater', font=('Arial',12,'bold'))
lbl.configure(bg='#ffffff')
lbl.pack() #grid(row=2,column=1)

lbl3 = Label(ventana,text = 'Herramienta de calculo y evaluación de intercambiadores de Calor.', font=('Arial',12))
lbl3.configure(bg='#ffffff')
lbl3.pack()
lbl4 = Label(ventana,text = 'Desarrollado por CENICAÑA todos sus derechos reservados.', font=('Arial',12))
lbl4.configure(bg='#ffffff')
lbl4.pack()

lbl2 = Label(ventana,text = 'Seleccione el tipo de Intercambiador con el cual desea trabajar:', font=('Arial',12))
lbl2.configure(bg='#ffffff')
lbl2.pack()

lbl2 = Label(ventana,text = 'Intercambiador de Calor - Coraza y Tubos', font=('Arial',12))
lbl2.configure(bg='#ffffff')
lbl2.place(x=10, y=100)
btn1 = Button(ventana, text = 'Coraza y Tubos Liquido/Liquido \n (Flujos -> Temp)', command = mifuncion, bg='#ffa600', font=('Arial',12))
btn1.place(x=25, y=280) #grid(row=3, column =0)
btn4 = Button(ventana, text = 'Coraza y Tubos Liquido/Liquido \n (Temp -> Flujos)', command = mifuncion4, bg='#ffa600', font=('Arial',12))
btn4.place(x=25, y=350) #grid(row=3, column =0)
img1 = PhotoImage(file='STHX.gif')
lbl_img1 = Label(ventana, image=img1, bd=0).place(x=30, y=130)
btn2 = Button(ventana, text = 'Intercambiador de Calor Placas', command = mifuncion2, bg='#ffa600', font=('Arial',12))
btn2.place(x=335, y=130) #grid(row=3, column =1)
img2 = PhotoImage(file='PHX1.gif')
lbl_img2 = Label(ventana, image=img2, bd=0).place(x=345, y=200)
btn3 = Button(ventana, text = 'Intercambiador de Calor con Aletas', command = mifuncion3, bg='#ffa600', font=('Arial',12))
btn3.place(x=625, y=130) # grid(row=3, column =2)
img3 = PhotoImage(file='AHX.gif')
lbl_img3 = Label(ventana, image=img3, bd=0).place(x=630, y=200)

lbl4 = Label(ventana,text = 'DDU', font=('Arial',8))
lbl4.configure(bg='#ffffff')
lbl4.place(x=860, y=880)

ventana.mainloop()
