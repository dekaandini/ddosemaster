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
from PyQt6.QtGui import QPixmap
from ui_mainddose import Ui_MainWindow # Sesuaikan nama class UI main Anda
import func1
import func2
import func3
import func4
# Import class UI popup Anda
from ui_invalpopup import Ui_Dialog1
from ui_gofpopup import Ui_Dialog2

#import func2

class gofpopup(QtWidgets.QDialog):
    def __init__(self, parent=None, data_fit=None, model_type="f1"):
        super(gofpopup, self).__init__(parent)
        self.ui = Ui_Dialog2()
        self.ui.setupUi(self)
        self.ui.buttonBox.accepted.connect(self.accept) 
        # Tombol Cancel akan mengembalikan QDialog.Rejected
        self.ui.buttonBox.rejected.connect(self.reject)
    # --- PENGISIAN FIELD PARAMETER (KOREKSI KEY) ---
        self.received_data = data_fit
        self.param_widgets ={
            "A1_opt": [self.ui.a1LineEdit, self.ui.a1Label],
            "cv_A1_fit": [self.ui.cVA1LineEdit, self.ui.cVA1Label],
            "A2_opt": [self.ui.a2LineEdit, self.ui.a2Label],
            "cv_A2_fit": [self.ui.cVA2LineEdit, self.ui.cVA2Label],
            "lambda1_opt": [self.ui.lambda1LineEdit, self.ui.lambda1Label],
            "cv_lambda1_fit": [self.ui.lambda2LineEdit, self.ui.lambda2Label],
            "lambda3_opt": [self.ui.lambda2LineEdit_2, self.ui.lambda2Label_2],
            "cv_lambda3_fit": [self.ui.cVLambda2LineEdit, self.ui.cVLambda2Label],
            "text_corrmat": [self.ui.textEdit_2,self.ui.label_17]
        }
        self.model_mapping = {
            "f1": ["A1_opt", "cv_A1_fit", "lambda1_opt", "cv_lambda1_fit","text_corrmat"],
            "f2": ["A1_opt", "cv_A1_fit", "lambda1_opt", "cv_lambda1_fit", "text_corrmat"],
            "f3": ["A1_opt", "cv_A1_fit", "A2_opt", "cv_A2_fit", "lambda1_opt", "cv_lambda1_fit", "text_corrmat"],
            "f4": ["A1_opt", "cv_A1_fit", "lambda1_opt", "cv_lambda1_fit", "lambda3_opt", "cv_lambda3_fit", "text_corrmat"]
        }
        if self.received_data:
            self.setup_display(model_type)

    def setup_display(self, model_type):
        active_params = self.model_mapping.get(model_type, [])

        # Sembunyikan semua terlebih dahulu, lalu tampilkan yang dibutuhkan
        for param, widgets in self.param_widgets.items():
            is_visible = param in active_params
            for w in widgets:
                w.setVisible(is_visible)
            
            # Jika terlihat, isi datanya
            if is_visible:
                val = self.received_data.get(param, 0)
                if isinstance(val, (int, float)):
        # Jika angka, gunakan format 3 angka di belakang koma
                    widgets[0].setText(f"{val:.3f}")
                else:
                    # Jika teks (seperti text_corrmat), tampilkan langsung tanpa format .3f
                    widgets[0].setText(str(val))
                # Menggunakan .get() agar tidak error jika key tidak ada
                

        """self.ui.a1LineEdit.setText(f"{self.received_data['A1_opt']:.3f}")
            self.ui.cVA1LineEdit.setText(f"{self.received_data['cv_A1_fit']:.3f}")
            self.ui.a2LineEdit.setText(f"{self.received_data['A2_opt']:.3f}")
            self.ui.cVA2LineEdit.setText(f"{self.received_data['cv_A2_fit']:.3f}")
            self.ui.lambda1LineEdit.setText(f"{self.received_data['lambda1_opt']:.3f}")
            self.ui.lambda2LineEdit.setText(f"{self.received_data['cv_lambda1_fit']:.3f}")
            self.ui.textEdit_2.setPlainText(self.received_data['text_corrmat'])
            #self.ui.textEdit_2.setText(self.received_data['text_corrmat']:.2f)"""
            #self.ui.a1LineEdit.setText(str(self.received_data["A1_opt"]))
            #f"{self.fit_r4:.3f}"
class invalpopup(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        # Load file UI dialog Anda
        self.ui = Ui_Dialog1()
        self.ui.setupUi(self)
        self.ui.buttonBox.accepted.connect(self.accept) 
        # Tombol Cancel akan mengembalikan QDialog.Rejected
        self.ui.buttonBox.rejected.connect(self.reject)
    
    def get_values(self):
        # Mengambil nilai dari QLineEdit di dialog
        return {
            "A1": float(self.ui.a1LineEdit.text()or 0),
            "A2": float(self.ui.a2LineEdit.text() or 0),
            "L1": float(self.ui.lambda1LineEdit.text() or 0),
            "L2": float(self.ui.lambda2LineEdit.text() or 0)
        }
    
    def set_values(self, data):
        # Menampilkan nilai hasil hitungan otomatis ke dialog
        self.ui.a1LineEdit.setText(str(round(data['A1'], 4)))
        self.ui.a2LineEdit.setText(str(round(data['A2'], 4)))
        self.ui.lambda1LineEdit.setText(str(round(data['L1'], 4)))
        self.ui.lambda2LineEdit.setText(str(round(data['L2'], 4)))
        # ... dst    
class MainDeka(QtWidgets.QMainWindow):
    # 1. Definisikan Signal Kustom (membawa dua list: x dan y)
 

    def __init__(self):
        super(MainDeka, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle("D-TIAC")
        self.t_data = None
        self.y_data = None
        self.fit_r1 = None # Inisialisasi juga agar tidak AttributeErr
        self.fitt_r1 = None
        self.fittt_r1 = None
        self.fitttt_r1 = None
        self.ui.cb_pg1.activated.connect(self.handle_param_guess1)
        self.ui.cb_pg1_2.activated.connect(self.handle_param_guess2)
        self.ui.cb_pg1_3.activated.connect(self.handle_param_guess3)
        self.ui.cb_pg1_4.activated.connect(self.handle_param_guess4)
        
        #FUNGSI TOMBOL
        #self.ui.fit1_button.clicked.connect(self.run_f1)
        
        # 3. Koneksikan Tombol Plot ke fungsi pengambilan data
 
        self.ui.gof1_button.clicked.connect(lambda: self.open_gof_dialog("f1"))
        self.ui.pushButton_8.clicked.connect(lambda: self.open_gof_dialog("f2"))
        self.ui.pushButton_12.clicked.connect(lambda: self.open_gof_dialog("f3"))
        self.ui.pushButton_15.clicked.connect(lambda: self.open_gof_dialog("f4"))

#TOMBOL FIG1
        self.ui.fit1_button.clicked.connect(self.fit_button1)
        self.ui.res1_button.clicked.connect(self.res_button1)
        #self.ui.gof1_button.clicked.connect(self.gof_button1)
        self.ui.cal1_button.clicked.connect(self.cal_button1)
#TOMBOL FIG2     
        self.ui.pushButton_7.clicked.connect(self.fit_button2)
        self.ui.pushButton_9.clicked.connect(self.res_button2)
        #self.ui.pushButton_12.clicked.connect(self.gof_button3)
        self.ui.pushButton_10.clicked.connect(self.cal_button2)

#TOMBOL FIG3     
        self.ui.pushButton_11.clicked.connect(self.fit_button3)
        self.ui.pushButton_13.clicked.connect(self.res_button3)
        #self.ui.pushButton_12.clicked.connect(self.gof_button3)
        self.ui.pushButton_17.clicked.connect(self.cal_button3)

#TOMBOL FIG4     
        self.ui.pushButton_14.clicked.connect(self.fit_button4)
        self.ui.pushButton_16.clicked.connect(self.res_button4)
        #self.ui.pushButton_12.clicked.connect(self.gof_button3)
        self.ui.pushButton_18.clicked.connect(self.cal_button4)

        """self.ui.gof1_button.clicked.connect(self.display_gof1)
        self.ui.pushButton_7.clicked.connect(self.jalankan_semua1)
        self.ui.pushButton_9.clicked.connect(self.res_plot2)
        self.ui.pushButton_10.clicked.connect(self.display_tiac_tia2)
        self.ui.pushButton_8.clicked.connect(self.display_gof2)"""
        self.ui.pushButton_5.clicked.connect(self.tampilkan_gambar)
        self.ui.pushButton_4.clicked.connect(self.load_excel_to_table)
        self.ui.pushButton_6.clicked.connect(self.conv_satuan)
        self.ui.pushButton.clicked.connect(self.calculate_waicc)
        self.ui.pushButton_3.clicked.connect(self.absorbed_dose)
    def handle_plot_button(self):
        # 4. Logika Mengambil Data dari Table Widget
        x_vals = []
        y_vals = []

        row_count = self.ui.input_table.rowCount() # Ganti nama tableWidget sesuai UI Anda

        for row in range(row_count):
            # Ambil item kolom 0 (Time) dan kolom 1 (% IA)
            item_x = self.ui.input_table.item(row, 0)
            item_y = self.ui.input_table.item(row, 1)

            # Pastikan cell tidak kosong
            if item_x and item_y and item_x.text() and item_y.text():
                try:
                    x = float(item_x.text())
                    y = float(item_y.text())
                    x_vals.append(x)
                    y_vals.append(y)
                except ValueError:
                    print(f"Data baris {row} bukan angka valid.")
        
        # 5. Emit (Kirim) Data jika tidak kosong
        if x_vals and y_vals:
            self.t_data = np.array(x_vals)
            self.y_data = np.array(y_vals)
            self.fit_r1 = func1.curve_fitting1(self.t_data, self.y_data)
            self.fitt_r1 = func3.curve_fitting3(self.t_data, self.y_data)
            self.fittt_r1 = func2.curve_fitting2(self.t_data, self.y_data)
            self.fitttt_r1 = func4.curve_fitting4(self.t_data, self.y_data)

        else:
            print("Data tabel kosong atau tidak valid.")

    def load_excel_to_table(self):
    # 1. Buka File Dialog untuk memilih file
        #options = QtWidgets.QFileDialog.Options()
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Pilih File Excel", "", "Excel Files (*.xlsx *.xls)"
        )

        if file_path:
            try:
            # 2. Baca file menggunakan pandas
                df = pd.read_excel(file_path)
                self.t_data = df.iloc[:, 0].values  # Mengambil semua baris kolom pertama
                self.y_data = df.iloc[:, 1].values  # Mengambil semua baris kolom kedua
            # Pastikan kolom sesuai (misal: kolom 0 adalah Time, kolom 1 adalah % IA)
            # 3. Iterasi data dan masukkan ke QTableWidget
            # self.tableWidget adalah nama variabel tabel Anda di GUI
                for row_index, row_data in df.iterrows():
                # Pastikan tidak melebihi jumlah baris di tabel GUI Anda
                    if row_index < self.ui.input_table.rowCount():
                    
                    # Kolom 1 (Time h) - index 0 di tabel
                        val_x = str(row_data.iloc[0]) 
                        self.ui.input_table.setItem(row_index, 0, QtWidgets.QTableWidgetItem(val_x))
                    
                    # Kolom 2 (% IA) - index 1 di tabel
                        val_y = str(row_data.iloc[1])
                        self.ui.input_table.setItem(row_index, 1, QtWidgets.QTableWidgetItem(val_y))
            
                print(f"Data berhasil dimuat!")
            except Exception as e:
                print(f"Terjadi kesalahan: {e}")
    #def run_f1(self):
    def open_gof_dialog(self, model_type):
    # Pilih data berdasarkan model_type
        data_map = {
            "f1": self.fit_r1,
            "f2": self.fittt_r1, # Sesuaikan nama variabel data Anda
            "f3": self.fitt_r1,
            "f4": self.fitttt_r1
        }
        
        current_data = data_map.get(model_type)
        
        if current_data is not None:
            dialog = gofpopup(self, data_fit=current_data, model_type=model_type)
            dialog.exec()
        else:
            print(f"Data untuk {model_type} belum tersedia. Lakukan fitting terlebih dahulu.")    
    def conv_satuan(self):
        try:
            input_ia = self.ui.IAedit.text()
            ia_cal = float(input_ia)
            ia_conv = ia_cal * 37
            self.ui.lineEdit.setText(f"{ia_conv:.1f}")
        except ValueError:
            self.ui.lineEdit.setText("Input Error!")

    def tampilkan_gambar(self):
        file_path = "Development(D-TIAC)/ginjalst_organ.png"
        pixmap = QPixmap(file_path)
        #lebar = self.ui.label.width()
        #tinggi = self.ui.label.height()
        lebar = 280
        tinggi = 170
        self.ui.label.setFixedSize(lebar, tinggi)

        pixmap_scaled = pixmap.scaled(
            lebar,
            tinggi,
            QtCore.Qt.AspectRatioMode.KeepAspectRatio,
            QtCore.Qt.TransformationMode.SmoothTransformation
        )
        self.ui.label.setPixmap(pixmap_scaled)
        self.ui.label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
    # Di dalam class MainDeka, bagian __init__
    

    def handle_param_guess1(self):
        pilihan = self.ui.cb_pg1.currentText()
        dialog = invalpopup()

        if pilihan == "Manual":
            # Mode Manual: User mengisi, lalu kita ambil nilainya
            if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                self.user_params = dialog.get_values()
             

        elif pilihan == "Automatic":
            # Mode Otomatis: Hitung dulu di func1, lalu tampilkan hasilnya di dialog
            # Contoh panggil fungsi hitung otomatis di func1
            auto_val = func1.calculate_auto1(self.t_data, self.y_data)
            
            # Tampilkan hasil hitungan ke dialog agar user bisa melihat
            dialog.set_values(auto_val)
            #dialog.setDisabled(True) # Opsional: buat read-only agar user tidak mengubah
            dialog.exec()
            
            self.user_params = auto_val
    #TOMBOL FIT
    def fit_button1(self):
        if self.ui.groupBox_2.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1.currentText()
                self.fit_r1 = func1.curve_fitting1(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )

                # Sekarang Anda bisa mengambil nilainya dengan aman
                self.fit_r2 = self.fit_r1["t_fit"]
                self.fit_r3 = self.fit_r1["y_fit"]
                self.fit_r4 = self.fit_r1["r_squared"]
                self.fit_r5 = self.fit_r1["aicc_aw"]
                
                self.ui.outr21.setText(f"{self.fit_r4:.3f}")
                self.ui.aICc1.setText(f"{self.fit_r5:.3f}")
                self.ui.aICc1LineEdit.setText(f"{self.fit_r5:.3f}")
                plt.scatter(self.t_data, self.y_data, label='Data Eksperimen')
                plt.plot(self.fit_r2, self.fit_r3, 'r-', label='Fitting Monoeksponensial')
                plt.xlabel('Waktu (t)')
                plt.ylabel('y')
                plt.legend()
                plt.grid(True)
                plt.show()
                            # ... lanjut ke plotting ...
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            self.ui.outr21.setText("")
            #self.curve_fitting()

    def res_button1(self):
        if self.ui.groupBox_2.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1.currentText()
                self.fit_r1 = func1.curve_fitting1(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.fit_r6 = self.fit_r1["residuals"]
                plt.figure()
                plt.scatter(self.t_data, self.fit_r6)
                plt.axhline(0, color='red', linestyle='--')
                plt.xlabel('Waktu (t)')
                plt.ylabel('Residuals')
                plt.title('Plot Residuals vs Waktu')
                plt.grid()
                plt.show()
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
             


    def gof_button1(self):
        if self.ui.groupBox_2.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1.currentText()
                self.fit_r1 = func1.curve_fitting1(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.popup = gofpopup(parent=self, data_fit=self.fit_r1)
                self.popup.show()
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
            #self.curve_fitting()
    

    def cal_button1(self):
        if self.ui.groupBox_2.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1.currentText()
                self.fit_r1 = func1.curve_fitting1(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.fit_r7 = self.fit_r1["tia_value"]
                self.fit_r8 = self.fit_r1["tia_unc"]
                self.fit_r9 = self.fit_r1["tiac_value"]
                self.fit_r10 = self.fit_r1["tiac_unc"]
                self.ui.tIAUncertaintyLineEdit.setText(f"{self.fit_r7:.3f} ± {self.fit_r8:.3f}")
                self.ui.tIACUncertaintyLineEdit.setText(f"{self.fit_r9:.3f} ± {self.fit_r10:.3f}")
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
            #self.curve_fitting()
        
    def handle_param_guess2(self):
        pilihan = self.ui.cb_pg1_2.currentText()
        dialog = invalpopup()

        if pilihan == "Manual":
            # Mode Manual: User mengisi, lalu kita ambil nilainya
            if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                self.user_params = dialog.get_values()
             

        elif pilihan == "Automatic":
            # Mode Otomatis: Hitung dulu di func1, lalu tampilkan hasilnya di dialog
            # Contoh panggil fungsi hitung otomatis di func1
            auto_val = func2.calculate_auto2(self.t_data, self.y_data)
            
            # Tampilkan hasil hitungan ke dialog agar user bisa melihat
            dialog.set_values(auto_val)
            #dialog.setDisabled(True) # Opsional: buat read-only agar user tidak mengubah
            dialog.exec()
            
            self.user_params = auto_val
    #TOMBOL FIT
    def fit_button2(self):
        if self.ui.groupBox_3.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_2.currentText()
                self.fittt_r1 = func2.curve_fitting2(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )

                # Sekarang Anda bisa mengambil nilainya dengan aman
                self.fit_r2 = self.fittt_r1["t_fit"]
                self.fit_r3 = self.fittt_r1["y_fit"]
                self.fit_r4 = self.fittt_r1["r_squared"]
                self.fit_r5 = self.fittt_r1["aicc_aw"]
                
                self.ui.outr21_2.setText(f"{self.fit_r4:.3f}")
                self.ui.aICcLineEdit_2.setText(f"{self.fit_r5:.3f}")
                self.ui.aICc2LineEdit.setText(f"{self.fit_r5:.3f}")
                plt.scatter(self.t_data, self.y_data, label='Data Eksperimen')
                plt.plot(self.fit_r2, self.fit_r3, 'r-', label='Fitting Bieksponensial-iv')
                plt.xlabel('Waktu (t)')
                plt.ylabel('y')
                plt.legend()
                plt.grid(True)
                plt.show()
                            # ... lanjut ke plotting ...
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            self.ui.outr21_2.setText("")
            #self.curve_fitting()

    def res_button2(self):
        if self.ui.groupBox_3.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_2.currentText()
                self.fittt_r1 = func2.curve_fitting2(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.fit_r6 = self.fittt_r1["residuals"]
                plt.figure()
                plt.scatter(self.t_data, self.fit_r6)
                plt.axhline(0, color='red', linestyle='--')
                plt.xlabel('Waktu (t)')
                plt.ylabel('Residuals')
                plt.title('Plot Residuals vs Waktu')
                plt.grid()
                plt.show()
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
             


    def gof_button2(self):
        if self.ui.groupBox_3.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_2.currentText()
                self.fittt_r1 = func2.curve_fitting2(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.popup = gofpopup(parent=self, data_fit=self.fittt_r1)
                self.popup.show()
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
            #self.curve_fitting()
    

    def cal_button2(self):
        if self.ui.groupBox_3.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_2.currentText()
                self.fittt_r1 = func2.curve_fitting2(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.fit_r7 = self.fittt_r1["tia_value"]
                self.fit_r8 = self.fittt_r1["tia_unc"]
                self.fit_r9 = self.fittt_r1["tiac_value"]
                self.fit_r10 = self.fittt_r1["tiac_unc"]
                self.ui.tIAUncertaintyLineEdit_2.setText(f"{self.fit_r7:.3f} ± {self.fit_r8:.3f}")
                self.ui.tIACUncertaintyLineEdit_2.setText(f"{self.fit_r9:.3f} ± {self.fit_r10:.3f}")
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
            #self.curve_fitting()
        

    def handle_param_guess3(self):
        pilihan = self.ui.cb_pg1_3.currentText()
        dialog = invalpopup()

        if pilihan == "Manual":
            # Mode Manual: User mengisi, lalu kita ambil nilainya
            if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                self.user_params = dialog.get_values()
             

        elif pilihan == "Automatic":
            # Mode Otomatis: Hitung dulu di func1, lalu tampilkan hasilnya di dialog
            # Contoh panggil fungsi hitung otomatis di func1
            auto_val = func3.calculate_auto3(self.t_data, self.y_data)
            
            # Tampilkan hasil hitungan ke dialog agar user bisa melihat
            dialog.set_values(auto_val)
            #dialog.setDisabled(True) # Opsional: buat read-only agar user tidak mengubah
            dialog.exec()
            
            self.user_params = auto_val
    #TOMBOL FIT
    def fit_button3(self):
        if self.ui.groupBox_4.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_3.currentText()
                self.fitt_r1 = func3.curve_fitting3(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )

                # Sekarang Anda bisa mengambil nilainya dengan aman
                self.fit_r2 = self.fitt_r1["t_fit"]
                self.fit_r3 = self.fitt_r1["y_fit"]
                self.fit_r4 = self.fitt_r1["r_squared"]
                self.fit_r5 = self.fitt_r1["aicc_aw"]
                
                self.ui.outr21_3.setText(f"{self.fit_r4:.3f}")
                self.ui.aICcLineEdit_3.setText(f"{self.fit_r5:.3f}")
                self.ui.aICc3LineEdit.setText(f"{self.fit_r5:.3f}")
                plt.scatter(self.t_data, self.y_data, label='Data Eksperimen')
                plt.plot(self.fit_r2, self.fit_r3, 'r-', label='Fitting Bieksponensial')
                plt.xlabel('Waktu (t)')
                plt.ylabel('y')
                plt.legend()
                plt.grid(True)
                plt.show()
                            # ... lanjut ke plotting ...
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            self.ui.outr21_3.setText("")
            #self.curve_fitting()

    def res_button3(self):
        if self.ui.groupBox_4.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_3.currentText()
                self.fitt_r1 = func3.curve_fitting3(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.fit_r6 = self.fitt_r1["residuals"]
                plt.figure()
                plt.scatter(self.t_data, self.fit_r6)
                plt.axhline(0, color='red', linestyle='--')
                plt.xlabel('Waktu (t)')
                plt.ylabel('Residuals')
                plt.title('Plot Residuals vs Waktu')
                plt.grid()
                plt.show()
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
             


    def gof_button3(self):
        if self.ui.groupBox_4.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_3.currentText()
                self.fitt_r1 = func3.curve_fitting3(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.popup = gofpopup(parent=self, data_fit=self.fitt_r1)
                self.popup.show()
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
            #self.curve_fitting()
    

    def cal_button3(self):
        if self.ui.groupBox_4.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_3.currentText()
                self.fitt_r1 = func3.curve_fitting3(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.fit_r7 = self.fitt_r1["tia_value"]
                self.fit_r8 = self.fitt_r1["tia_unc"]
                self.fit_r9 = self.fitt_r1["tiac_value"]
                self.fit_r10 = self.fitt_r1["tiac_unc"]
                self.ui.tIAUncertaintyLineEdit_3.setText(f"{self.fit_r7:.3f} ± {self.fit_r8:.3f}")
                self.ui.tIACUncertaintyLineEdit_3.setText(f"{self.fit_r9:.3f} ± {self.fit_r10:.3f}")
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
            #self.curve_fitting()
        
    def handle_param_guess4(self):
        pilihan = self.ui.cb_pg1_4.currentText()
        dialog = invalpopup()

        if pilihan == "Manual":
            # Mode Manual: User mengisi, lalu kita ambil nilainya
            if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                self.user_params = dialog.get_values()
             

        elif pilihan == "Automatic":
            # Mode Otomatis: Hitung dulu di func1, lalu tampilkan hasilnya di dialog
            # Contoh panggil fungsi hitung otomatis di func1
            auto_val = func4.calculate_auto4(self.t_data, self.y_data)
            
            # Tampilkan hasil hitungan ke dialog agar user bisa melihat
            dialog.set_values(auto_val)
            #dialog.setDisabled(True) # Opsional: buat read-only agar user tidak mengubah
            dialog.exec()
            
            self.user_params = auto_val
    #TOMBOL FIT
    def fit_button4(self):
        if self.ui.groupBox_5.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_4.currentText()
                self.fitttt_r1 = func4.curve_fitting4(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )

                # Sekarang Anda bisa mengambil nilainya dengan aman
                self.fit_r2 = self.fitttt_r1["t_fit"]
                self.fit_r3 = self.fitttt_r1["y_fit"]
                self.fit_r4 = self.fitttt_r1["r_squared"]
                self.fit_r5 = self.fitttt_r1["aicc_aw"]
                
                self.ui.outr21_4.setText(f"{self.fit_r4:.3f}")
                self.ui.aICcLineEdit_4.setText(f"{self.fit_r5:.3f}")
                self.ui.aICc4LineEdit.setText(f"{self.fit_r5:.3f}")
                plt.scatter(self.t_data, self.y_data, label='Data Eksperimen')
                plt.plot(self.fit_r2, self.fit_r3, 'r-', label='Fitting Monoeksponensial')
                plt.xlabel('Waktu (t)')
                plt.ylabel('y')
                plt.legend()
                plt.grid(True)
                plt.show()
                            # ... lanjut ke plotting ...
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            self.ui.outr21_4.setText("")
            #self.curve_fitting()

    def res_button4(self):
        if self.ui.groupBox_5.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_4.currentText()
                self.fitttt_r1 = func4.curve_fitting4(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.fit_r6 = self.fitttt_r1["residuals"]
                plt.figure()
                plt.scatter(self.t_data, self.fit_r6)
                plt.axhline(0, color='red', linestyle='--')
                plt.xlabel('Waktu (t)')
                plt.ylabel('Residuals')
                plt.title('Plot Residuals vs Waktu')
                plt.grid()
                plt.show()
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
             


    def gof_button4(self):
        if self.ui.groupBox_5.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_4.currentText()
                self.fitttt_r1 = func4.curve_fitting4(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.popup = gofpopup(parent=self, data_fit=self.fitttt_r1)
                self.popup.show()
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
            #self.curve_fitting()
    

    def cal_button4(self):
        if self.ui.groupBox_5.isChecked():
            if self.t_data is not None and self.y_data is not None:
            # Panggil fungsi dari func1.py
                error_type = self.ui.errorModelDataComboBox.currentText()
                objective_func = self.ui.objectiveFunctionComboBox.currentText()
                initial_value = self.ui.cb_pg1_4.currentText()
                self.fit_tttr1 = func4.curve_fitting4(self.t_data, self.y_data, err_mod = error_type, obj_func = objective_func, in_val = self.user_params )
                self.fit_r7 = self.fitttt_r1["tia_value"]
                self.fit_r8 = self.fitttt_r1["tia_unc"]
                self.fit_r9 = self.fitttt_r1["tiac_value"]
                self.fit_r10 = self.fitttt_r1["tiac_unc"]
                self.ui.tIAUncertaintyLineEdit_4.setText(f"{self.fit_r7:.3f} ± {self.fit_r8:.3f}")
                self.ui.tIACUncertaintyLineEdit_4.setText(f"{self.fit_r9:.3f} ± {self.fit_r10:.3f}")
            else:
                print("Error: Silakan klik Plot atau input data terlebih dahulu!")
        else:
            return
            #self.curve_fitting()
        

    def calculate_waicc (self):
        aicc1 = self.fit_r1["aicc_aw"]
        tiac1 = self.fit_r1["tiac_value"]
        tiac_unc1 = self.fit_r1["tiac_unc"]
        aicc2 = self.fittt_r1["aicc_aw"]
        tiac2 = self.fittt_r1["tiac_value"]
        tiac_unc2 = self.fittt_r1["tiac_unc"]
        aicc3 = self.fitt_r1["aicc_aw"]
        tiac3 = self.fitt_r1["tiac_value"]
        tiac_unc3 = self.fitt_r1["tiac_unc"]
        aicc4 = self.fitttt_r1["aicc_aw"]
        tiac4 = self.fitttt_r1["tiac_value"]
        tiac_unc4 = self.fitttt_r1["tiac_unc"]
        
        data={
            'Model': ['f1', 'f2', 'f3', 'f4'],
            'AICc': [aicc1, aicc2, aicc3, aicc4],
            'TIAC': [tiac1, tiac2, tiac3, tiac4],
            'SE_TIAC': [tiac_unc1, tiac_unc2, tiac_unc3, tiac_unc4]            
        }
        df = pd.DataFrame(data)

        aicc_min = df['AICc'].min()
        df['Delta_i'] = df['AICc'] - aicc_min
        df['Likelihood'] = np.exp(-df['Delta_i'] / 2)
        total_likelihood = df['Likelihood'].sum()

        df['w_AICc'] = df['Likelihood'] / total_likelihood
        self.avg_TIAC = (df['w_AICc'] * df['TIAC']).sum()

        #Menghitung SE(TIAC)
        variance_term = (df['TIAC'] - self.avg_TIAC)**2
        se_squared_term = df['SE_TIAC']**2

        sqrt_term = np.sqrt(se_squared_term + variance_term)
        df['se_component'] = df['w_AICc'] * sqrt_term

        avg_SE_TIAC = df['se_component'].sum()

        waicc1 = df.loc[0, 'w_AICc']
        waicc2 = df.loc[1, 'w_AICc']
        waicc3 = df.loc[2, 'w_AICc']
        waicc4 = df.loc[3, 'w_AICc']
        
        self.ui.wAICc1LineEdit.setText(f"{waicc1:.2f}")
        self.ui.wAICc2LineEdit.setText(f"{waicc2:.2f}")
        self.ui.wAICc3LineEdit.setText(f"{waicc3:.2f}")
        self.ui.wAICc4LineEdit.setText(f"{waicc4:.2f}")
        self.ui.tIACAvgLineEdit.setText(f"{self.avg_TIAC:.4f}")
        self.ui.standarErrorSELineEdit.setText(f"{avg_SE_TIAC:.6f}")

    def absorbed_dose(self):
        tiac_final = self.avg_TIAC * 3600
        s_value_gamma = 1.3444533e-6
        s_value_beta = 6.5988141e-5
        s_tot = s_value_gamma +  s_value_beta
        absorbeddose = tiac_final * s_tot
        self.ui.abdose.setText(f"{absorbeddose:.6f}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainDeka()
    window.show()
    sys.exit(app.exec())
            #self.send_data_signal.emit(x_vals, y_vals)
            
            # Jika receiver adalah popup window dan belum muncul, munculkan sekarang
            
            #self.receiver_widget.show() 
