import numpy as np
import matplotlib.pyplot as pltcd
import pandas as pd
import statistics
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
from scipy.optimize import minimize
from scipy import integrate
from numdifftools import Hessian
from scipy.stats import linregress
from PyQt6 import QtWidgets, QtCore
from PyQt6.QtCore import Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt6.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QLabel # Import class receiver tadi
#FUNGSI MONOEXP
#from ui_gofpopup import Ui_Dialog

t_phys = 159.463  #t1/2 177 Lu-DOTATATE (jam)
lambda2 = np.log(2)/t_phys #lambda2 = lambda physics = koefisien peluruhan fisik (jam^-1)
#initial_guess = [in_val["A1"], in_val["L1"]]
#initial_guess = [1.0, 0.02]
bounds = [(0, None), (0, None)]  # A1 dan lambda1 >= 0


def calculate_auto1 (t_data, y_data):
    log_y_data = np.log(y_data)
    slope, intercept = np.polyfit(t_data, log_y_data, 1)
    A1_initial = np.exp(intercept)
    lambda1_initial = -slope - lambda2
    initial_guess = [A1_initial, lambda1_initial]
    return {
        "A1": A1_initial,
        "A2": 0, # Sesuaikan jika modelnya f1 (single exponential)
        "L1": lambda1_initial,
        "L2": 0
    }

def mono_exponential(t, A1, lambda1) :
    return A1*np.exp(-(lambda1 + lambda2) *t)


def objective_function1 (params, t_data, y_data, err_mod, obj_func) :
    #A, B, C = 0, 0.01, 2
    #sigma_data = np.sqrt(A + B * y_data**C)
    
    #params, params_covariance = curve_fit(bi_exponentialiv, t_data, y_data, p0=initial_guess) #pcov dievaluasi terhadap y data
    A1, lambda1 = params
    #params = params
    #A1_mpop = 2.232 #diambil dari data Pak Barron dengan mengecualikan data parameter pasien 1
    #lambda1_mpop = 0.038 #diambil dari data Pak Barron dengan mengecualikan data parameter pasien 1
    n = len(y_data)
    phi = 3.14
    #Ada if disini untuk opsi error model data
    if err_mod == "FSD 10%":
        A, B, C = 0, 0.01, 2
        sigma_data = np.sqrt(A + B * y_data**C)
    else :
        np.random.seed(42)
        jumlah_data = len(y_data)
        bobot_awal = np.random.rand(jumlah_data)
        sigma_data = bobot_awal / np.sum(bobot_awal)
    
    A1_mpop = 2.232 #diambil dari data Pak Barron dengan mengecualikan data parameter pasien 1
    lambda1_mpop = 0.038 
    y_pred = mono_exponential(t_data, A1, lambda1)
    residuals = (y_data-y_pred)/sigma_data
    residualss = residuals**2
    ofs3 = np.sum(residualss)
    v = ofs3 / n
    vhi = 2*phi*(sigma_data**2)
    ofs1 = np.log(vhi)
    ofs2 = np.log(v)

    if obj_func == "SSE":
        oftot = ofs1 + ofs2 + ofs3
        return np.sum(oftot)
    else:
     #diambil dari data Pak Barron dengan mengecualikan data parameter pasien 1
        #tambahan komponen Bayessian (wj samakan dulu dengan sigma_data)
        wjbayes = 2*phi*(sigma_data**2)
        ofs4 = np.log(wjbayes)
        bayes1 = (A1 - A1_mpop)
        pjbayes1 = (bayes1 / sigma_data)**2
        bayes2 = (lambda1 - lambda1_mpop)
        pjbayes2 = (bayes2 / sigma_data)**2
        ofbay = ofs4 + pjbayes1 + pjbayes2
        oftot = ofs1 + ofs2 + ofs3 + ofbay  
        #ada if sebelum return untuk opsi OF
        return np.sum(oftot)        
def curve_fitting1 (t_data, y_data, err_mod, obj_func, in_val) :
    #t_phys = 159.463  #t1/2 177 Lu-DOTATATE (jam)
    #lambda2 = np.log(2)/t_phys #lambda2 = lambda physics = koefisien peluruhan fisik (jam^-1)

    #ada if disini untuk opsi initial value parameter
    
    #Ada if disini untuk opsi error model data
    #A, B, C = 0, 0.01, 2
    #sigma_data = np.sqrt(A + B * y_data**C)
    initial_guess = [in_val["A1"], in_val["L1"]]
    params,params_covariance = curve_fit(mono_exponential, t_data, y_data, p0=initial_guess) #pcov dievaluasi terhadap y data
    result = minimize(
        fun= objective_function1,
        x0= initial_guess,
        args=(t_data, y_data, err_mod, obj_func),
        method='trust-constr',
        bounds= bounds,
        options={'verbose': 1}  # Cetak proses iterasi
    )
    A1_opt, lambda1_opt = result.x 
    t_fit = np.linspace(0, 180, 100)
    y_fit = mono_exponential(t_fit, A1_opt, lambda1_opt)
    res = y_data - mono_exponential(t_data, *params)  
    ss_res = np.sum(res**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (ss_res / ss_tot)
    #ui.outr21_2.setText(f"{r_squared:.3f}")
    #plt.plot(t_fit, y_fit, 'r-', label='Fitting Monoeksponensial')
    #plt.show()
    #HITUNG AICC
    n = len(y_data)  # Jumlah data points
    k = 2  
    obj_func = objective_function1(result.x, t_data, y_data, err_mod, obj_func)
    #ada if disini untuk opsi absolute/realtive weighting dengan cek nilai N dan K
    aicc_aw = obj_func - 2*k + (2*k*(k+1) / (n-k-1)) #absolute
    #PLOT RESIDUALS
    residuals = y_data - mono_exponential(t_data, *params)

    #HITUNG GOF RESULTS
    variance = np.diag(params_covariance)
    std_error = np.sqrt(variance)
    A1_mean = params[0]
    lambda1_mean = params[1]
    cv_A1_fit = (std_error[0] / A1_mean) * 100
    cv_lambda1_fit = (std_error[1] / lambda1_mean) * 100
    outer_std_devs = np.outer(std_error, std_error)
    corrmat = params_covariance / outer_std_devs 
    fix_corrmat = np.round(corrmat, 2)
    text_corrmat = str(fix_corrmat)

    #HITUNG TIA DAN TIAC
    suk1 = A1_opt / (lambda1_opt + lambda2)
    H = Hessian(lambda p: objective_function1(p, t_data, y_data, err_mod, obj_func))(result.x)
    cov_matrix = np.linalg.inv(H)
    grad_tia = np.array([
    1 / (lambda1_opt + lambda2),  # dTIAC/dA1
    -A1_opt / (lambda1_opt + lambda2)**2 # dTIAC/dlambda1
    ])
    grad_tiac = np.array([
    1 / (lambda1_opt + lambda2)/100,  # dTIAC/dA1
    -A1_opt / (lambda1_opt + lambda2)**2/100 # dTIAC/dlambda1
    ])
    tia_unc = np.sqrt(grad_tia.T @ cov_matrix @ grad_tia)
    tiac_unc = np.sqrt(grad_tiac.T @ cov_matrix @ grad_tiac)
    tia_value = suk1
    tiac_value = tia_value /100

    return {
        "t_fit": t_fit,
        "y_fit": y_fit,
        "r_squared": r_squared,
        "aicc_aw": aicc_aw,
        "residuals": residuals,
        "tia_value": tia_value,
        "tia_unc": tia_unc,
        "tiac_value": tiac_value,
        "tiac_unc": tiac_unc,
        "A1_opt": A1_opt,
        "cv_A1_fit": cv_A1_fit,
        "lambda1_opt": lambda1_opt,
        "cv_lambda1_fit": cv_lambda1_fit,
        "text_corrmat": text_corrmat

    }

    #ui.aICcLineEdit_2.setText(f"{aicc_aw:.2f}")
    #ui.aICc2LineEdit.setText(f"{aicc_aw:.2f}")   


    """plt.figure()
    plt.scatter(t_data, residuals)
    plt.axhline(0, color='red', linestyle='--')
    plt.xlabel('Waktu (t)')
    plt.ylabel('Residuals')
    plt.title('Plot Residuals vs Waktu')
    plt.grid()
    plt.show() """
"""def calculate_tiac1() :
    t_phys = 159.463  #t1/2 177 Lu-DOTATATE (jam)
    lambda2 = np.log(2)/t_phys #lambda2 = lambda physics = koefisien peluruhan fisik (jam^-1)
    A1_opt, lambda1_opt = result.x 
    suk1 = 1 / (lambda1_opt + lambda2)
    H = Hessian(lambda p: objective_function1(p, t_data, y_data, sigma_data))(result.x)
    cov_matrix = np.linalg.inv(H)
    grad_tia = np.array([
    1 / (lambda1_opt + lambda2),  # dTIAC/dA1
    -A1_opt / (lambda1_opt + lambda2)**2 # dTIAC/dlambda1
    ])
    grad_tiac = np.array([
    1 / (lambda1_opt + lambda2)/100,  # dTIAC/dA1
    -A1_opt / (lambda1_opt + lambda2)**2/100 # dTIAC/dlambda1
    ])
    tia_unc = np.sqrt(grad_tia.T @ cov_matrix @ grad_tia)
    tiac_unc = np.sqrt(grad_tiac.T @ cov_matrix @ grad_tiac)
    tia_value = suk1
    tiac_value = tia_value /100
    return{
        "tia_value": tia_value,
        "tia_unc": tia_unc,
        "tiac_value": tiac_value,
        "tiac_unc": tiac_unc
    }
    #ui.tIAUncertaintyLineEdit.setText(f"{tia_value:.2f} ± {tia_unc:.2f}")
    #ui.tIACUncertaintyLineEdit.setText(f"{tiac_value:.2f} ± {tiac_unc:.2f}")
    #A1_opt, lambda1_opt = result.x """

"""def hitung_qcpar() :
    residuals = y_data - bi_exponentialiv(t_data, *params)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (ss_res / ss_tot)
    n = len(y_data)  # Jumlah data points
    k = 2
    obj_func = objective_function(params, t_data, y_data, sigma_data)
    aicc_aw = obj_func - 2*k + (2*k*(k+1) / (n-k-1)) 
    outr21.setText(f"{r_squared:.3f}") # Format dengan 3 desimal
    aiCc1.setText(f"{aicc_aw:.2f}")"""

"""def calculate_gof_data1(params, params_covariance):
    #gof_dialog = gofpopup()
    A1_opt, lambda1_opt = result.x 
    A1, lambda1 = params
    params_covariance = params_covariance
    variance = np.diag(params_covariance)
    std_error = np.sqrt(variance)
    A1_mean = params[0]
    lambda1_mean = params[1]
    cv_A1_fit = (std_error[0] / A1_mean) * 100
    cv_lambda1_fit = (std_error[1] / lambda1_mean) * 100
    outer_std_devs = np.outer(std_error, std_error)
    corrmat = params_covariance / outer_std_devs 
    fix_corrmat = np.round(corrmat, 1)
    text_corrmat = str(fix_corrmat)
    return {
        "A1_opt": A1_opt,
        "cv_A1_fit": cv_A1_fit,
        "lambda1_opt": lambda1_opt,
        "cv_lambda1_fit": cv_lambda1_fit,
        "text_corrmat": text_corrmat
    }
    gof_dialog.ui.a1LineEdit.setText(f"{A1_opt:.2f}%")
    gof_dialog.ui.cVA1LineEdit.setText(f"{cv_A1_fit:.2f}%")
    gof_dialog.ui.lambda1LineEdit.setText(f"{lambda1_opt:.2f}%")
    gof_dialog.ui.lambda2LineEdit.setText(f"{cv_lambda1_fit:.2f}%")
    gof_dialog.ui.textEdit.setText(text_corrmat)
    gof_dialog.exec()  # Tampilkan dialog secara modal"""

   #TOMBOL FIT
"""def jalankan_semua1():
    # Di sini kita panggil keduanya berurutan
    #handle_plot_button()

    curve_fitting1()
    minimize_of1()
    
def display_tiac_tia1():
    calculate_tiac1()
    
    # TIA dan TIAC (Asumsi sudah berupa string dengan Uncertainty)
    #lineEdit_TIA.setText(tia_value)
    #lineEdit_TIAC.setText(tiac_value)
def display_gof1():
    calculate_gof_data1(params,params_covariance)""" 
