from PyQt5 import QtGui, QtCore  # (the example applies equally well to PySide)
import pyqtgraph as pg 
from fortprof import profile as prof
from fortmag import magnetosphere as mag
import readmodel as rm
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.constants import year,pi,G
from scipy.optimize import fsolve
# from os import remove as remove_file
import os
import glob
import numpy as np

def get_T(grid, nh, ne, T_start, Rstar):
		log_T_for_fit = np.array([3, 3.7, 3.8, 3.9, 4, 4.2, 4.6, 4.9, 5.4])
		log_lamb_for_fit = np.array([-44.4, -28.3, -26.0, -24.5, -23.6, -22.6, -21.8, -21.2, -21.2])
		
		log_lamb_log_T = InterpolatedUnivariateSpline(log_T_for_fit, log_lamb_for_fit, k=3)
		def log_lamb(T):
				return log_lamb_log_T(np.log(T)/np.log(10))
		def lamb(T):
				return 10**log_lamb(T)

		# M_sun = 1.989e33

		# T_plot = np.linspace(5100, 10000, num = 100, endpoint = True)
		# lamb_plot = log_lamb(T_plot)
		# plt.plot(T_plot, lamb_plot, 'k')
		# plt.show()
		# quit()

		R_sun = 6.957e10
		
		R_star = Rstar*R_sun
		# M_star = Mstar*M_sun
		
		rms = grid[:,:,0]
		tets = grid[:,:,1]
		gridr,gridtet = nh.shape
		distance = np.array([ rms[i]*np.sin(tets[i,:])**2 for i in range(len(rms)) ])*R_star
		# print(distance[-1,2])
		# if(heat == 0):
		heat_coef = lamb(T_start)*nh[0,0]*ne[0,0]*distance[0,0]**3
		# else:
				# heat_coef = heat
		def solve_T_bal(grid,nh,ne):

				def find_T(T_start):
						heat_rate = heat_coef/distance**3
						# print('f=',f)

						def func(T, heat, nh,ne):
								# print(lamb(T))
								return lamb(T)*nh*ne-heat

						T = np.zeros(nh.shape)
						for i in range(gridr):
								for j in range(gridtet):
										# print(i,j, heat_rate[i,j], nh[i,j], ne[i,j])
										T[i,j]= fsolve(func, T_start, (heat_rate[i,j], nh[i,j], ne[i,j]))
										# print(T[i,j], lamb(T[i,j]))
						# print(T)
						return T


				T = find_T(T_start)

				return T

		ind_max = 0
		new_max = 1
		T = np.zeros(nh.shape)
		T = solve_T_bal(grid,nh,ne)
		while ind_max != new_max:
				ind_max = new_max
				Tm = T[-1,:]
				T_max = max(Tm)
				new_max = np.argmax(Tm)
				nh_max = nh[-1, ind_max]
				ne_max = ne[-1, ind_max]
				dis_max = distance[-1, ind_max]
				# print(T_max, nh_max, ne_max, dis_max)
				T_max = T_start
				heat_coef = lamb(T_max)*nh_max*ne_max*dis_max**3
				# print(heat_coef)
				T = solve_T_bal(grid,nh,ne)
				pass
		# T = solve_T_bal(grid,nh,ne)
		# Tm = max(T[-1,:])
		# print(new_max, np.argmax(T[-1,:]))
		return T

model_dir = 'models/'
populations_dir = 'models/popul/'
model_data_dir = 'models/data/'
hotspot_data_dir = 'models/hotspot/'
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


class Star:
	def __init__(self, mass, radius, temperature, name = '', vrot = 0):
		self.mass = mass
		self.radius = radius
		self.name = name
		self.temperature = temperature
		self.vrot = vrot

class Field:
	def __init__(self, mdot, geom_type, first_border, second_border, cut_in, cut_out):
		self.type = geom_type
		self.mdot = mdot
		self.first_border = first_border
		self.second_border = second_border
		self.cut_in = cut_in
		self.cut_out = cut_out

class Model:
	def __init__(self, name = '', new = False, data = ''):
		self.name = name
		self.star = None
		self.field = None
		self.populations = name
		self.populations_calculated = False
		self.populations_deattached = False
		self.hotspot_temp = []
		self.hotspot_width = []
		self.hotspot_pos = []
		if not new:
			self._fromfile()
		else:
			self._new(data)
		print(self.hotspot_file_name)
		self.init_hotspot_file()
		self.initial_populations = (self.populations == self.name)
		self.populations_file = populations_dir + self.populations + '_popul.dat'
		self.check_populations()
		if self.populations_calculated:
			nh, ne, levels = self.get_populations()
		else:
			if self.field.type == 'dipole':
				self.grid_accuracy = tuple([5, 20])
			if self.field.type == 'cone':
				self.grid_accuracy = tuple([20])
			self.grid_calculated = False
			self.temperature_calculated = False
		self.grid_step = np.array([[0, 0.5], [1, 0.5]])

	def _fromfile(self):
		try:
			population_parameters = rm.read_parameters(self.name)
		except rm.NoSuchModelError as err:
			print("Can't find data for model "+err.model_name)
			population_parameters = None
		Rstar = population_parameters['Rstar']
		Mstar = population_parameters['Mstar']
		Tstar = population_parameters['Tstar']
		star_name = self.name.split('_')[0]
		self.star = Star(Mstar, Rstar, Tstar, name = star_name, vrot = 0)
		first_border = population_parameters['first_border']
		second_border = population_parameters['second_border']
		in_cut = population_parameters['in_cut']
		out_cut = population_parameters['out_cut']
		Mdot = population_parameters['Mdot']
		field_type = population_parameters['field_type']
		self.field = Field(Mdot, field_type, first_border, second_border, in_cut, out_cut)
		self.populations = population_parameters['populations']
		hotspot = population_parameters['hot_spot']
		self.Tmax = population_parameters['Tmax']
		self.additional_temperature = np.zeros(self.grid_accuracy)
		# if hotspot:
		#   pos = 0.0
		#   width = population_parameters['Dhot']
		#   temp = population_parameters['Thot']
		#   self.add_hotspot(pos, temp, width)

	def _new(self, data):
		population_parameters = rm.read_parameters_from_string(self.name, data)
		Rstar = population_parameters['Rstar']
		Mstar = population_parameters['Mstar']
		Tstar = population_parameters['Tstar']
		star_name = self.name.split('_')[0]
		self.star = Star(Mstar, Rstar, Tstar, name = star_name, vrot = 0)
		first_border = population_parameters['first_border']
		second_border = population_parameters['second_border']
		in_cut = population_parameters['in_cut']
		out_cut = population_parameters['out_cut']
		Mdot = population_parameters['Mdot']
		field_type = population_parameters['field_type']
		self.field = Field(Mdot, field_type, first_border, second_border, in_cut, out_cut)
		self.populations = population_parameters['populations']
		self.hotspot_file_name = population_parameters['hotspot']
		# print('model init: ', hotspot_file_name)
		self.Tmax = population_parameters['Tmax']
		# self.additional_temperature = np.zeros(self.grid_accuracy)

	def init_hotspot_file(self):
		print(self.hotspot_file_name)
		if self.hotspot_file_name != '':
			self.hotspot_file = hotspot_data_dir + self.hotspot_file_name + '_hotspot.dat'
			self.add_hotspots()
		else:
			self.hotspot_file = ''

	def add_hotspots(self):
		# self.delete_hotspots()
		parameters = {'width' : self.hotspot_width, 'temp' : self.hotspot_temp, 'pos' : self.hotspot_pos}
		print('add')
		if self.field.type == 'dipole':
			try:
				# print('add', self.hotspot_file)
				hot_file = open(self.hotspot_file, 'r')
			except FileNotFoundError:
				return
				print('adderr')
			nval = 0
			# print(hot_file.readlines())
			for line in hot_file.readlines():
				# print(line)
				if line != '':
					line=line.replace(" ", '').replace("\n",'')
					values_str, parameter = tuple(line.split('#')[:2])
					values_words = values_str.split(',')
					print(values_words, parameter)
					if nval == 0:
						nval = len(values_words)
					elif nval != len(values_words):
						self.delete_hotspots()
						return
					values = [float(value_str) for value_str in values_words]
					parameters[parameter].extend(values)
			hot_file.close()

	def delete_hotspots(self):
		self.hotspot_pos = []
		self.hotspot_width = []
		self.hotspot_temp = []

	def init_temperature(self):
		if self.grid_calculated:
			Mdot = self.field.mdot
			Mstar = self.star.mass
			Rstar = self.star.radius
			Tstar = self.star.temperature
			nh, ne, ni = mag.calc_populations(self.grid, self.temperature, Mdot, Mstar, Rstar, Tstar, 0.0, 7, loc=True)
			self.temperature = get_T(self.grid, nh, ne, self.Tmax, Rstar)
			self.temperature_calculated = True

	def init_hotspots(self):
		if self.field.type == 'dipole':
			if self.grid_calculated:
				self.additional_temperature = np.zeros(self.grid_accuracy)
				gridr, gridtet = self.grid_accuracy
				d_grid = mag.calc_d_grid(gridtet, self.grid_step.transpose())
				for i, rm in enumerate(self.grid[:,0,0]):
					cos_theta = (1-1/rm**2)**0.5
					rm_d_grid = rm/2*(np.log(cos_theta+(1+cos_theta**2)**0.5)+cos_theta*(1+cos_theta**2)**0.5)
					# rm_d_grid = d_grid*rm_d_grid_constant
					for pos, width, temp in zip(self.hotspot_pos, self.hotspot_width, self.hotspot_temp):
						self.additional_temperature[i,:] = self.additional_temperature[i,:] + temp*np.exp(-np.abs(d_grid - pos)**2*rm_d_grid**2/width**2)

	def set_grid_accuracy(self, acc):
		self.grid_accuracy = acc
		self.grid_calculated = False

	def init_grid(self):
		if self.field.type == 'dipole':
			# print(self.grid_accuracy)
			gridr, gridtet = self.grid_accuracy
			rmi, rmo = self.field.first_border, self.field.second_border
			print('init grid step: ', self.grid_step)
			d_grid = mag.calc_d_grid(gridtet, self.grid_step.transpose())
			# print(d_grid)
			if np.sum(np.abs(d_grid)/gridtet) < 1e-6: 
				d_grid = np.linspace(0, 1, gridtet)
			self.grid = mag.init_grid(gridr, rmi, rmo, mtet=gridtet, grid_d = d_grid)
			self.grid_calculated = True
			self.temperature = np.full(self.grid_accuracy, self.Tmax)
			self.additional_temperature = np.zeros(self.grid_accuracy)
			self.temperature_calculated = False
			# self.init_temperature()
			# print(self.grid[0,:,1])

	def get_populations(self):
		grid_accuracy_found = False
		if self.populations_calculated:
			if self.field.type == 'dipole':
				irm = 0; itet = 0
				popul_file = open(self.populations_file, 'r')
				for line in popul_file.readlines():
					line = line.replace('\n', '')
					line = line.strip()
					# print(line.strip().replace(' ', '_'))
					if(line != '' and line[0] == '#'):
						try:
							grid_accuracy_str = line[1:].split()
							self.grid_accuracy = tuple([int(acc_str) for acc_str in grid_accuracy_str])
							grid_accuracy_found = True
							self.additional_temperature = np.zeros(self.grid_accuracy)
							temperature = np.zeros(self.grid_accuracy)
							nh = np.zeros(self.grid_accuracy); ne = np.zeros(self.grid_accuracy)
							mr, mtet = self.grid_accuracy
							levels = np.zeros((mr, mtet, 7))
							self.grid = np.zeros((mr, mtet, 2))
						except:
							pass
					elif grid_accuracy_found:
						if line != '':
							line_data_str = line.split()
							line_data = [float(data_str) for data_str in line_data_str]
							# print(line_data)
							self.grid[irm, itet, :] = [line_data[0], line_data[1]]
							temperature[irm, itet] = line_data[2]
							nh[irm, itet] = line_data[3]
							ne[irm, itet] = line_data[-1]
							levels[irm, itet, :] = line_data[4:11]
							itet += 1
						else:
							irm += 1
							itet = 0
				if grid_accuracy_found:
					self.grid_calculated = True
					self.temperature_calculated = True
					self.temperature = temperature
					# self.hotspot_file_name = ''
					self.init_hotspot_file()
					self.delete_hotspots()
				else:
					self.grid_calculated = False
					self.temperature_calculated = False
				return nh, ne, levels
			else:
				return [], [], []

	def check_populations(self):
		self.populations_calculated = True
		try:
			f = open(self.populations_file, 'r')
			f.close()
		except FileNotFoundError:
			self.populations_calculated = False

	def calculate(self):
		if not self.populations_calculated:
			print(self.name) 
			self.init_grid()
			self.init_temperature()
			self.delete_hotspots()
			self.add_hotspots()
			self.init_hotspots()
			temp = self.temperature + self.additional_temperature
			Mdot = self.field.mdot
			Mstar = self.star.mass
			Rstar = self.star.radius
			Tstar = self.star.temperature
			# if not self.grid_calculated:
			nh, ne, ni = mag.calc_populations(self.grid, temp, Mdot, Mstar, Rstar, Tstar, 0.0, 
																									15, model='temp', loc=False)
			self.populations_calculated = True
			os.rename('models/temp_popul.dat', self.populations_file)
		




	# def set_parameters:
	#   self.delete_hotspots()
	#   population_parameters = rm.read_parameters_from_string(self.name, data)
	#   Rstar = population_parameters['Rstar']
	#   Mstar = population_parameters['Mstar']
	#   Tstar = population_parameters['Tstar']
	#   star_name = self.name.split('_')[0]
	#   self.star = Star(Mstar, Rstar, Tstar, name = star_name, vrot = 0)
	#   first_border = population_parameters['first_border']
	#   second_border = population_parameters['second_border']
	#   in_cut = population_parameters['in_cut']
	#   out_cut = population_parameters['out_cut']
	#   Mdot = population_parameters['Mdot']
	#   field_type = population_parameters['field_type']
	#   self.field = Field(Mdot, field_type, first_border, second_border, in_cut, out_cut)
	#   self.populations = population_parameters['populations']
	#   self.hotspot_file_name = population_parameters['hotspot']
	#   # print('model init: ', hotspot_file_name)
	#   self.Tmax = population_parameters['Tmax']
	#   # self.additional_temperature = np.zeros(self.grid_accuracy)
	#   if self.hotspot_file_name != '':
	#     self.hotspot_file = hotspot_data_dir + self.hotspot_file_name + '_hotspot.dat'
	#     self.add_hotspots()
	#   else:
	#     self.hotspot_file = ''


class TextBox(QtGui.QWidget):
	def __init__(self, label, submit_label = 'OK'):
		QtGui.QWidget.__init__(self)
		self.label = QtGui.QLabel(label)
		self.text_line = QtGui.QLineEdit()
		self.submit_btn = QtGui.QPushButton(submit_label)
		layout = QtGui.QGridLayout()
		layout.setMargin(0)
		self.setLayout(layout)
		layout.addWidget(self.label, 0,0)
		layout.addWidget(self.text_line, 0,1)
		layout.addWidget(self.submit_btn, 0,2)

class EditModel(QtGui.QWidget):
	def __init__(self, parent = None):
		super(EditModel, self).__init__()
		self.parent = parent
		self.model_index = -1
		layout = QtGui.QGridLayout()
		layout.setMargin(0)
		# layout.setSpacing(0)
		layout.setRowStretch(0,0)
		layout.setRowStretch(1,1)
		self.setLayout(layout)
		self.model_name = ''
		self.model = None
		self.model_name_textbox = TextBox('Model', submit_label = 'New')
		self.model_name_textbox.text_line.setMinimumWidth(100)
		self.model_plain_text = QtGui.QTextEdit()
		# self.new_model_btn = QtGui.QPushButton('New')
		self.save_btn = QtGui.QPushButton('Add/Save')
		self.status_line = QtGui.QLabel()
		self.status_line.setAlignment(QtCore.Qt.AlignLeft)
		layout.addWidget(self.model_name_textbox, 0, 0, 1, 4)
		# layout.addWidget(self.new_model_btn, 0, 3)
		layout.addWidget(self.model_plain_text, 1, 0, 3, 4)
		layout.addWidget(self.status_line, 4, 0, 1, 3)
		layout.addWidget(self.save_btn, 4, 3)
		self.model_name_textbox.submit_btn.clicked.connect(self.create_model_data)
		# self.new_model_btn.clicked.connect(self.create_model_data)
		self.save_btn.clicked.connect(self.save_model)
 
	def open_model_data(self, forced = False):
		model_name = self.model_name_textbox.text_line.text()
		if (not forced) and self.model_name == model_name:
			self.set_status_line("Model "+model_name+" already opened", content = 'err')
			return
		self.model_name = model_name
		if self.model_name == '':
			self.set_status_line("Model name is empty, please set", content = 'err')
			return
		parameters_file_name = self.model_name + '_data.dat'
		try:
				parameters_file = open(model_data_dir+parameters_file_name, 'r')
				model_text = parameters_file.read()
				self.model_plain_text.setText(model_text)
		except FileNotFoundError:
				self.set_status_line("Can't find file "+parameters_file_name, content = 'err')
				self.model_name = ''
				return
		self.set_status_line('Opened model '+model_name, content = 'msg')
		self.model_index = -1
		try:
			model_data = self.model_plain_text.toPlainText()
			self.model = Model(self.model_name, new = True, data = model_data)
			uneditable = self.model.initial_populations and self.model.populations_calculated
			print(uneditable)
			self.model_plain_text.setReadOnly(uneditable)
		except rm.ImportantParameterNotSet as err:
			pass
		except rm.InvalidParameterValue as err:
			pass

		parameters_file.close()


	def create_model_data(self):
		model_name = self.model_name_textbox.text_line.text()
		if self.model_name == model_name:
			self.set_status_line("Please change model name first", content = 'err')
			return
		self.model_name = model_name
		if self.model_name == '':
			self.set_status_line("Model name is empty, please set", content = 'err')
			return
		try:
				parameters_file_name = self.model_name + '_data.dat'
				parameters_file = open(model_data_dir+parameters_file_name, 'r')
				self.set_status_line("Model "+model_name+" already exist", content = 'err')
				return
		except FileNotFoundError:
				pass       
		# parameters_file = open('models/'+parameters_file_name, '2')
		model_text = self.model_plain_text.toPlainText()
		self.model_index = -1
		if model_text == '':
			model_text = ('# Parameters for model '+self.model_name+'\n'+
					'float_value # Mstar\n'+ 'float_value # Tstar\n'+ 'float_value # Rstar\n'+
					'float_value # Mdot\n'+ 'float_value # Tmag\n'+ 'string_value (dipole or cone) # field type\n'+
					'# borders\n'+'float_value # first\n'+'float_value # second\n')
			self.model_plain_text.setText(model_text)
			self.set_status_line('Created template for model ' + self.model_name, content = 'msg')
		else:
			model_text_lines = model_text.split('\n')
			model_text_lines[0] = '# Parameters for model '+self.model_name
			model_text = '\n'.join(model_text_lines)
			self.model_plain_text.setText(model_text)
			self.set_status_line('Changed model name to ' + self.model_name, content = 'msg')
		uneditable = False
		self.model_plain_text.setReadOnly(uneditable)
		print('model_name:', self.model_name)

	def save_model(self):
		print('model_name:', self.model_name)
		if self.model_name == '':
			self.set_status_line("Model name is empty, please set", content = 'err')
			return
		model_data = self.model_plain_text.toPlainText()
		print('model_name:', self.model_name)
		try:
			self.model = Model(self.model_name, new = True, data = model_data)
		except rm.ImportantParameterNotSet as err:
			self.set_status_line('Parameter '+err.parameter+' is important. Use "# '+err.parameter_key+'" to set',
																content = 'err')
			return
		except rm.InvalidParameterValue as err:
			self.set_status_line('Invalid value for parameter ' + err.parameter_key, content = 'err')
			return
		print('model_name:', self.model_name)
		if self.model_index == -1:
			self.model_index = self.parent.add_model_to_list(self.model)
		elif not(self.model.initial_populations and self.model.populations_calculated):
			del_model = Main.models[self.model_index]
			Main.models[self.model_index] = self.model
			del del_model
			self.parent.update_list()
		print('model_name:', self.model_name)
		self.parent.set_chosen_model(self.model_index)
		parameters_file_name = self.model_name + '_data.dat'
		parameters_file = open(model_data_dir+parameters_file_name, 'w')
		parameters_file.write(model_data)
		self.set_status_line('Saved to ' + parameters_file_name, content = 'msg')
		self.parent.update_model_files()
		self.parent.load_model_view(self.model_index)
		# self.parent.edit_model_by_index(self.model_index)

	def load_model(self, model_index, forced = False):
		self.model_index = model_index
		self.model = Main.models[model_index]
		# self.model_name = self.model.name
		self.model_name_textbox.text_line.setText(self.model.name)
		parameters_file_name = self.model_name + '_data.dat'
		self.open_model_data(forced = forced)
		self.model_index = model_index

	def set_status_line(self, message, content = 'msg'):
		contents_style = {'msg' : "color : #000000",
											'err' : "color : #aa0000"}
		self.status_line.setText(message)
		self.status_line.setStyleSheet(contents_style[content])

class ModelManager(QtGui.QWidget):
	def __init__(self, parent):
		super(ModelManager, self).__init__()
		layout = QtGui.QGridLayout()
		layout.setMargin(0)
		self.setLayout(layout)
		self.parent = parent
		self.delete_model_btn = QtGui.QPushButton('Delete model')
		self.delete_file_btn = QtGui.QPushButton('Delete file')
		self.delete_model_btn.clicked.connect(self.delete_model)
		self.delete_file_btn.clicked.connect(self.delete_model_file)
		self.model_list = QtGui.QListWidget()
		# self.model_list.setIconSize(QtCore.QSize(10,10))
		self.model_files_list = QtGui.QListWidget()
		self.model_files_status_line = QtGui.QLabel('')
		self.model_files_status_line.setAlignment(QtCore.Qt.AlignLeft)
		self.model_list_header = QtGui.QLabel('Loaded models')
		self.model_files_list_header = QtGui.QLabel('Model files')
		self.calculate_btn = QtGui.QPushButton('Calculate')
		self.calculate_btn.clicked.connect(self.calculate_models)
		self.model_list.itemClicked.connect(self.edit_model)
		self.model_list.itemDoubleClicked.connect(self.model_populations_attachment)
		# self.model_list.itemSelectionChanged.connect(self.model_selected)
		# self.model_files_list.itemDoubleClicked.connect(self.add_model_from_file)
		self.model_files_list.itemClicked.connect(self.load_model_from_file)

		layout.addWidget(self.delete_file_btn, 0, 1)
		layout.addWidget(self.model_files_list_header, 0, 0)
		layout.addWidget(self.model_files_list, 1, 0, 3, 2)
		layout.addWidget(self.delete_model_btn, 0, 3)
		layout.addWidget(self.calculate_btn, 0, 2, 1, 1)
		layout.addWidget(self.model_list, 1, 2, 4, 2)
		layout.addWidget(self.model_files_status_line, 4, 0, 2, 1)
		layout.setRowStretch(3,1)
		layout.setRowStretch(4,0)
		layout.setRowMinimumHeight(4,25)
		self.edit_model_widget = EditModel(self)
		layout.addWidget(self.edit_model_widget, 0, 4, 5, 3)
		self.update_model_files()
	
	def calculate_models(self):
		models = Main.models
		for index, model in enumerate(models):
			self.set_chosen_model(index)
			model.calculate()
		self.update_list()


	def model_populations_attachment(self, item):
		index = self.model_list.selectedIndexes()[0].row()
		model = Main.models[index]
		if model.populations_calculated:
			model.populations_calculated = False
			model.populations_deattached = True
			model.grid_calculated = False
			model.temperature_calculated = False
		elif model.populations_deattached:
			model.populations_calculated = True
			nh,ne,levels = model.get_populations()
		self.update_list()
		self.edit_model_by_index(index)


	def edit_model(self, item):
		index = self.model_list.selectedIndexes()[0].row()
		self.edit_model_by_index(index)

	def edit_model_by_index(self, index):
		self.edit_model_widget.load_model(index, forced = True)
		self.load_model_view(index)

	def load_model_view(self, index):
		self.parent.model_view.load_model(index)

	def delete_model_file(self):
		index = self.model_files_list.selectedIndexes()[0].row()
		# item = self.model_files_list.selectedItems()[0]
		item = self.model_files_list.item(index)
		model_name = item.text()
		model_names = [model.name for model in Main.models]
		if not (model_name in model_names):
			self.model_files_list.takeItem(index)
			model_file = model_data_dir+model_name+'_data.dat'
			os.remove(model_file)
			self.set_files_status_line('Deleted '+model_name+' file', content='msg')
		else:
			self.set_files_status_line('Remove '+model_name+' from models first', content='err')

	def set_files_status_line(self, message, content = 'msg'):
		contents_style = {'msg' : "color : #000000",
											'err' : "color : #aa0000"}
		self.model_files_status_line.setText(message)
		self.model_files_status_line.setStyleSheet(contents_style[content])

	def delete_model(self):
		selected = self.model_list.selectedIndexes()
		if len(selected) < 1: return
		index = self.model_list.selectedIndexes()[0].row()
		del Main.models[index]
		self.update_list()
		self.edit_model_widget.model_index = -1

	def update_model_files(self):
		self.model_files_list.clear()
		model_files_paths = glob.glob(model_data_dir + '*_data.dat')
		model_files = [model_file_path.split('/')[-1] for model_file_path in model_files_paths]
		model_names = ['_'.join(model_file.split('_')[:-1]) for model_file in model_files]
		for model_name in model_names:
			self.model_files_list.addItem(model_name)

	def load_model_from_file(self, item):
		model_name = item.text()
		# self.edit_model_widget.model_name = model_name
		self.edit_model_widget.model_name_textbox.text_line.setText(model_name)
		self.edit_model_widget.open_model_data()

	def add_model_from_file(self, item):
		model_name = item.text()
		self.load_model_from_file(item)
		# self.edit_model(item)
		self.edit_model_widget.save_model() 
		self.set_files_status_line('Added '+model_name+' to models', content='msg')


	def update_list(self):
		self.model_list.clear()
		for model in Main.models:
			item = QtGui.QListWidgetItem(model.name, self.model_list)
			if model.populations_calculated:
				icon = QtGui.QIcon('icons/yes.png')
			else:
				if model.populations_deattached:
					icon = QtGui.QIcon('icons/yesbutactuallyno.png')
				else:
					icon = QtGui.QIcon('icons/no.png')
			item.setIcon(icon)
		self.parent.model_view.reset()
			# item.setIconSize
			# print(icon.isNull())
	
	def set_chosen_model(self, index):
		print('set chosen ', index)
		self.model_list.setCurrentRow(index)
		# self.edit_model_by_index(index)

	def add_model_to_list(self, model):
		model_names = [model.name for model in Main.models]
		if not (model.name in model_names):
			index = len(Main.models)
			Main.models.append(model)
			item = QtGui.QListWidgetItem(model.name, self.model_list)
			if model.populations_calculated:
				icon = QtGui.QIcon('icons/yes.png')
			else:
				icon = QtGui.QIcon('icons/no.png')
			# icon = QtGui.QIcon.fromTheme("edit-undo")
			item.setIcon(icon)
		else:
			index = -1
		return index
		# icon.pixmap(512).save('test.png')
		# self.model_list.addItem(model.name)

class ModelGridView(QtGui.QWidget):
	def __init__(self, parent):
		self.model = None
		self.parent = parent
		super(ModelGridView, self).__init__()
		self.grid_accuracy_edit = QtGui.QLineEdit()
		self.grid_accuracy_edit.setValidator(QtGui.QRegExpValidator(
												QtCore.QRegExp('^\d{1,2},\s*\d{1,2}')))
		self.grid_accuracy_label = QtGui.QLabel('Grid accuracy')
		self.save_btn = QtGui.QPushButton('Save grid')
		self.load_btn = QtGui.QPushButton('Load grid')
		d_grid_plot = pg.GraphicsWindow()
		grid_plot = pg.GraphicsWindow()
		self.d_grid_plot = pg.PlotItem()
		self.d_grid_chosen = pg.ScatterPlotItem()
		self.d_grid_plot.disableAutoRange()
		self.d_grid_plot.setXRange(-0.1, 1.1)
		self.d_grid_plot.setXRange(0.0, 1.0)
		self.d_grid_plot.setMouseEnabled(x = False, y = False)
		self.d_grid_plot.setMenuEnabled(False)
		label_font = QtGui.QFont()
		label_font.setWeight(QtGui.QFont.Normal)
		self.dvl = pg.InfiniteLine(angle = 90, pen = pg.mkPen(color = 'r'), label = 'label')
		self.dvl.label.setPosition(0.1); self.dvl.label.setColor('r'); self.dvl.label.setFont(label_font)
		self.dhl = pg.InfiniteLine(angle = 0, pen = pg.mkPen(color = 'r'), label = 'label')
		self.dhl.label.setPosition(0.1); self.dhl.label.setColor('r'); self.dhl.label.setFont(label_font)
		self.grid_plot = pg.PlotItem()
		self.grid_plot.setMenuEnabled(False)
		self.grid_plot.setMouseEnabled(x = False, y = False)

		d_grid_plot.addItem(self.d_grid_plot)
		# self.d_grid_plot.addItem(self.dhl)
		# self.d_grid_plot.addItem(self.dvl)
		self.d_grid_plot.addItem(self.d_grid_chosen)
		grid_plot.addItem(self.grid_plot)
		self.saved_grid_step = None
		self.saved_acc = None
		self.is_grid_saved = False
		
		# self.grid_plot.scene().sigMouseClicked.connect(self.test)
		layout = QtGui.QGridLayout()
		self.setLayout(layout)
		layout.addWidget(self.grid_accuracy_edit, 0, 1)
		layout.addWidget(self.grid_accuracy_label, 0, 0)
		layout.addWidget(self.save_btn, 0, 2)
		layout.addWidget(self.load_btn, 0, 3)
		layout.addWidget(d_grid_plot, 1, 2, 1, 2)
		layout.addWidget(grid_plot, 1, 0, 1, 2)
		
		self.model_index = -1
		self.chosen_index = -1
		# self.model.grid_step = np.array([[0, 0.5], [1, 0.5]])
		self.d_grid_plot_line = pg.PlotDataItem(pen = pg.mkPen(0.0, width = 2))
		self.d_grid_plot.addItem(self.d_grid_plot_line)

	def wheelEvent(self, event):
		delta = event.angleDelta().y()/120*0.001
		model = Main.models[self.model_index]
		if self.chosen_index >= 0:
			model.grid_step[self.chosen_index,1] += delta
			self.d_grid_plot_line.setData(model.grid_step.transpose()[0,:], model.grid_step.transpose()[1,:])
			self.update_chosen()
			model.grid_calculated = False
			self.init_grid()

	def init_grid(self):
		model = Main.models[self.model_index]
		plot = self.grid_plot
		plot.plot([0], [0], clear = True)
		if self.model_index > -1:
			if not model.grid_calculated:
				model.init_grid()
			mr, mt = model.grid_accuracy
			grid_points_rm = model.grid[:,:,0].reshape((mt*mr))
			grid_points_tet = model.grid[:,:,1].reshape((mt*mr))
			grid_points_x = grid_points_rm*np.sin(grid_points_tet)**3
			grid_points_y = grid_points_rm*np.sin(grid_points_tet)**2*np.cos(grid_points_tet)
			plot.addItem(pg.ScatterPlotItem(grid_points_x, grid_points_y))
		# self.parent.grid_updated = True

	def update_chosen(self):
		if self.chosen_index >= 0:
			model = Main.models[self.model_index]
			self.d_grid_chosen.setData([model.grid_step[self.chosen_index, 0]], [model.grid_step[self.chosen_index, 1]])
			# model.grid_calculated = False
		else:
			self.d_grid_chosen.setData([], [])

	def d_grid_mousemove(self, event):
		# print(event)
		model = Main.models[self.model_index]
		pos = event
		plot = self.d_grid_plot
		vb = plot.vb
		x = [model.grid_step[0,0], model.grid_step[-1,0]]+list(model.grid_step[1:-1,0])
		y = [model.grid_step[0,1], model.grid_step[-1,1]]+list(model.grid_step[1:-1,1])
		if plot.sceneBoundingRect().contains(pos):
			mousePoint = vb.mapSceneToView(pos)
			xm = mousePoint.x(); ym = mousePoint.y()
			# print(x,y)
			self.dvl.setValue(xm); self.dvl.label.setText('x = {:3.2f}'.format(xm))
			self.dhl.setValue(ym); self.dhl.label.setText('y = {:3.2f}'.format(ym))
			if(ym < 0.5):
				self.dvl.label.setPosition(0.9)
			else:
				self.dvl.label.setPosition(0.1)
			if(xm < 0.5):
				self.dhl.label.setPosition(0.9)
			else:
				self.dhl.label.setPosition(0.1)
			if xm < 0: 
				min_ind = 0
				self.chosen_index = 0
			elif xm <= 1:
				if len(x) > 2:
					min_ind = np.argmin(abs(np.array(x)-xm)[2:])+2
					self.chosen_index = min_ind - 1
				else:
					self.chosen_index = -1
			else:
				min_ind = 1
				self.chosen_index = len(x)-1
			# xo = x[min_ind]; yo = y[min_ind]
			self.update_chosen()



		self.d_grid_plot.addItem(self.dhl)
		self.d_grid_plot.addItem(self.dvl)

	def d_grid_edit(self, event):
		model = Main.models[self.model_index]
		# print(event.pos(), event.)
		pos = event.scenePos()

		plot = self.d_grid_plot
		vb = plot.vb
		# print(self.model.grid_step)
		x = [model.grid_step[0,0], model.grid_step[-1,0]]+list(model.grid_step[1:-1,0])
		y = [model.grid_step[0,1], model.grid_step[-1,1]]+list(model.grid_step[1:-1,1])
		# print(plot.viewGeometry())
		if plot.sceneBoundingRect().contains(pos):
			model.grid_calculated = False
			mousePoint = vb.mapSceneToView(pos)
			x_m = mousePoint.x()
			y_m = mousePoint.y()
			# print(event.button())
			min_ind = 0
			if 0 < x_m < 1:
				if len(x) > 2:
					min_ind = np.argmin(abs(np.array(x)-x_m)[2:])+2
			if x_m <= 0:
				min_ind = 0
			if x_m > 1:
				min_ind = 1
			if event.button() == 1:
				if 0 < x_m < 1:
					x.append(x_m); y.append(y_m)
				if x_m <= 0:
					y[0] = y_m
				if x_m > 1:
					y[1] = y_m
			if event.button() == 2 and len(x) > 2 and min_ind != 0 and min_ind != 1:
				del x[min_ind]
				del y[min_ind]
				self.chosen_index = -1
			if event.button() == 4:
				y[min_ind] = y_m
			model.grid_step = np.array(list(zip(x,y)))
			model.grid_step = model.grid_step[model.grid_step[:,0].argsort()]
			# print(self.model.grid_step)
			# np.append(self.model.grid_step, np.array([[x],[y]]), axis = 1)
			# print(self.model.grid_step[0][:])
			# print(self.model.grid_step[0][0])
			# print(self.model.grid_step[:,0])
		self.d_grid_plot_line.setData(model.grid_step.transpose()[0,:], model.grid_step.transpose()[1,:])
		self.d_grid_mousemove(pos)
		self.init_grid()

	def save_grid(self):
		model = Main.models[self.model_index]
		self.saved_grid_step = model.grid_step
		self.saved_acc = model.grid_accuracy
		self.is_grid_saved = True
	
	def load_grid(self):
		model = Main.models[self.model_index]
		if self.is_grid_saved:
			print('loaded saved grid')
			model.grid_accuracy = self.saved_acc
			model.grid_step = self.saved_grid_step
			model.grid_calculated = False
			m_rm, m_thet = model.grid_accuracy

			self.grid_accuracy_edit.setText(str(m_rm) + ', ' + str(m_thet))
			# self
			self.d_grid_plot_line.setData(model.grid_step.transpose()[0,:], model.grid_step.transpose()[1,:])
			# self.d_grid_plot.addItem(self.d_grid_plot_line)
			# self.grid_btn.clicked.connect()
			self.init_grid()


	def load_model(self, model_index):
		self.grid_plot.plot([0], [0], clear = True)
		self.chosen_index = -1
		self.update_chosen()
		# if self.model_index > -1:
			# del_model = Main.models[self.model_index]
			# Main.models[self.model_index] = self.model
			# del del_model
		try:
			self.d_grid_plot.scene().sigMouseClicked.disconnect()
		except:
			pass
		try:
			self.d_grid_plot.scene().sigMouseMoved.disconnect()
		except:
			pass
		try:
			self.save_btn.clicked.disconnect()
		except:
			pass
		try:
			self.load_btn.clicked.disconnect()
		except:
			pass
		try:
			self.grid_accuracy_edit.editingFinished.disconnect()
		except:
			pass
		if model_index >= 0:
			self.model_index = model_index
			model = Main.models[self.model_index]
			if not model.populations_calculated:
				self.d_grid_plot.scene().sigMouseClicked.connect(self.d_grid_edit)
				self.d_grid_plot.scene().sigMouseMoved.connect(self.d_grid_mousemove)
				self.grid_accuracy_edit.editingFinished.connect(self.update_accuracy)
				self.save_btn.clicked.connect(self.save_grid)
				self.load_btn.clicked.connect(self.load_grid)
				

			field_type = model.field.type
			
			m_rm, m_thet = model.grid_accuracy

			self.grid_accuracy_edit.setText(str(m_rm) + ', ' + str(m_thet))
			# self
			self.d_grid_plot_line.setData(model.grid_step.transpose()[0,:], model.grid_step.transpose()[1,:])
			# self.d_grid_plot.addItem(self.d_grid_plot_line)
			# self.grid_btn.clicked.connect()
			self.init_grid()
		else: 

			self.grid_accuracy_edit.clear()
			# self.set_grid_status_line("Field type "+field_type+" is not supported yet", content = 'err')

	def set_grid_status_line(self, message, content = 'msg'):
		contents_style = {'msg' : "color : #000000",
											'err' : "color : #aa0000"}
		self.grid_accuracy_status.setText(message)
		self.grid_accuracy_status.setStyleSheet(contents_style[content])

	def update_accuracy(self):
		model = Main.models[self.model_index]
		accuracy_str = self.grid_accuracy_edit.text().replace(" ", '')
		accuracy_list = []
		for acc in accuracy_str.split(','):
			accuracy_list.append(int(acc))
		if(model != None):
			model.set_grid_accuracy(tuple(accuracy_list))
			self.init_grid()

class TemperatureWidget(QtGui.QWidget):
	def __init__(self, parent):
		self.parent = parent
		super(TemperatureWidget, self).__init__()
		layout = QtGui.QGridLayout()
		self.setLayout(layout)
		# self.model = None
		self.model_index = -1
		temp_plot = pg.GraphicsWindow()
		self.temp_plot = pg.PlotItem()
		self.temp_plot.setMenuEnabled(False)
		self.temp_plot.setMouseEnabled(x = False, y = False)
		self.temp_plot.setAutoVisible(y = True)
		self.temp_plot.enableAutoRange()
		self.temp_line = pg.PlotDataItem(pen = pg.mkPen(color = 'k', width = 2))
		self.add_temp_line = pg.PlotDataItem(pen = pg.mkPen(color = 'r', style=QtCore.Qt.DashLine))
		self.pcolormesh_temp = pg.widgets.MatplotlibWidget.MatplotlibWidget()
		fig = self.pcolormesh_temp.getFigure()
		win = fig.canvas.window() 
		toolbar = win.findChild(QtGui.QToolBar)
		toolbar.setVisible(False)
		self.pcolormesh_temp_plot = fig.subplots()
		self.pcolormesh_temp_plot.axis('equal')
		self.pcolormesh_temp_plot.set_adjustable('box')
		self.temp_plot.addItem(self.add_temp_line)
		self.temp_plot.addItem(self.temp_line)
		temp_plot.addItem(self.temp_plot)
		layout.addWidget(temp_plot, 0, 1, 5, 1)
		layout.addWidget(self.pcolormesh_temp, 0, 0, 5, 1)
		# self.name_text_line = QtGui.QLineEdit()
		self.pos_text_line = QtGui.QLineEdit()
		self.pos_text_line.setValidator(QtGui.QRegExpValidator(
												QtCore.QRegExp('^((0(\.\d+){,1})|(1(\.0+)){,1})(,\s*((0(\.\d+){,1})|(1(\.0+)){,1}))*')))
		self.width_text_line = QtGui.QLineEdit()
		self.width_text_line.setValidator(QtGui.QRegExpValidator(
												QtCore.QRegExp('^\d(\.\d+){,1}(,\s*(\d+(\.\d+){,1}))*')))
		self.temp_text_line = QtGui.QLineEdit()
		self.temp_text_line.setValidator(QtGui.QRegExpValidator(
												QtCore.QRegExp('^\d+(\.\d+){,1}(,\s*(\d+(\.\d+){,1}))*')))
		self.set_btn = QtGui.QPushButton('Set')
		# name_label = QtGui.QLabel('File name')
		pos_label = QtGui.QLabel('Position')
		width_label = QtGui.QLabel('Width')
		temp_label = QtGui.QLabel('Temperature')
		header_label = QtGui.QLabel('Hotspots parameters')
		layout.addWidget(header_label, 0, 2); layout.addWidget(self.set_btn, 0, 2)
		layout.addWidget(pos_label, 1, 2); layout.addWidget(self.pos_text_line, 1, 2)
		layout.addWidget(width_label, 2, 2); layout.addWidget(self.width_text_line, 2, 2)
		layout.addWidget(temp_label, 3, 2); layout.addWidget(self.temp_text_line, 3, 2)
		# layout.addWidget(name_label, 4, 1); layout.addWidget(self.name_text_line, 4, 2)
		layout.setColumnStretch(2, 0)
		layout.setColumnStretch(1, 2); layout.setColumnStretch(0, 4)
		# layout.setRowStretch(0, 0); layout.setRowStretch(1, 0); layout.setRowStretch(2, 0); layout.setRowStretch(3, 0)
		layout.setRowStretch(4, 10)

	def set_hotspots(self):
		# model = Main.models[self.model_index]
		# model.delete_hotspots()
		# model.set_hotspot_file()
		# if self.name_text_line.text() == '':
		#   return
		self.init_additional()

	def load_text_lines(self):
		parameters = {'width' : self.width_text_line, 
									'temp' : self.temp_text_line,
									'pos' : self.pos_text_line}
		model = Main.models[self.model_index]
		hotspot_file_name = model.hotspot_file
		try:
			hot_file = open(hotspot_file_name, 'r')
		except FileNotFoundError:
			return
		# self.name_text_line.setText(model.hotspot_file_name)
		for line in hot_file.readlines():
				# print(line)
				if line != '':
					line=line.replace(" ", '').replace("\n",'')
					values_str, parameter = tuple(line.split('#')[:2])
					parameters[parameter].setText(values_str)

	def load_model(self, model_index):
		self.pos_text_line.clear()
		self.width_text_line.clear()
		self.temp_text_line.clear()
		# self.name_text_line.clear()
		try:
			self.set_btn.clicked.disconnect()
		except:
			pass
		if model_index >= 0 and Main.models[model_index].grid_calculated:
			self.model_index = model_index
			self.load_text_lines()
			model = Main.models[self.model_index]
			print(model.name, model.hotspot_file)
			if not model.temperature_calculated:
				model.init_temperature()    
			self.update_temperature()
			model.init_hotspots()
			self.update_additional()
			if not model.populations_calculated:
				self.set_btn.clicked.connect(self.set_hotspots)
		else:
			self.temp_line.setData([],[])
			self.add_temp_line.setData([],[])
			fig = self.pcolormesh_temp.getFigure()
			fig.clear()
			self.pcolormesh_temp_plot = fig.subplots()
			self.pcolormesh_temp_plot.axis('equal')
			self.pcolormesh_temp_plot.set_adjustable('box')
		

	def init_additional(self):
		model = Main.models[self.model_index]
		print(model.name, model.hotspot_file)
		model.delete_hotspots()
		# file_name = self.pos_text_line.text()
		# model.set_hotspot_file(file_name)
		pos_text = self.pos_text_line.text()
		width_text = self.width_text_line.text()
		temp_text = self.temp_text_line.text()
		hotspot_file_name = model.hotspot_file
		print('hotspot file: ', hotspot_file_name)
		hotspot_file = open(hotspot_file_name, 'w')
		hotspot_file.write(pos_text + ' # pos\n' + 
											 width_text + ' # width\n' +
											 temp_text + ' # temp\n')
		hotspot_file.close()
		model.add_hotspots()
		model.init_hotspots()
		self.update_additional()

	def update_additional(self):
		def get_c(z):
			mx = len(z[:,0])
			my = len(z[0,:])
			c = np.zeros((mx,my))
			for i in range(mx-1):
				for j in range(my-1):
					c[i,j] = z[i,j] + z[i+1,j] + z[i,j+1] + z[i+1,j+1]       
			c = c/4
			return(c)
		model = Main.models[self.model_index]
		add_temp = model.additional_temperature
		print(add_temp[0,0], model.name, model.additional_temperature[0,0])
		plot = self.temp_plot
		mr, mt = model.grid_accuracy
		grid_points_rm = model.grid[:,:,0].reshape((mr*mt))
		grid_points_tet = model.grid[:,:,1].reshape((mr*mt))
		grid_points_r = grid_points_rm*np.sin(grid_points_tet)**2
		grid_x = model.grid[:,:,0]*np.sin(model.grid[:,:,1])**3
		grid_y = model.grid[:,:,0]*np.sin(model.grid[:,:,1])**2*np.cos(model.grid[:,:,1])
		temp = model.temperature + add_temp
		temp_points = temp.reshape((mr*mt))
		fig = self.pcolormesh_temp.getFigure()
		fig.clear()
		self.pcolormesh_temp_plot = fig.subplots()
		self.pcolormesh_temp_plot.axis('equal')
		self.pcolormesh_temp_plot.set_adjustable('box')
		temp_map = self.pcolormesh_temp_plot.pcolormesh(grid_x,grid_y,get_c(temp), cmap = 'Greys')
		cbar = fig.colorbar(temp_map, ax = self.pcolormesh_temp_plot, format='%.1f')
		self.pcolormesh_temp.draw()
		plot.setYRange(0.9*np.amin(temp_points), 1.1*np.amax(temp_points))
		connect = np.full((mr*mt), 1); connect[mt-1::mt] = 0
		self.add_temp_line.setData(grid_points_r, temp_points, connect = connect)

	def update_temperature(self):
		model = Main.models[self.model_index]
		plot = self.temp_plot
		mr, mt = model.grid_accuracy
		grid_points_rm = model.grid[:,:,0].reshape((mr*mt))
		grid_points_tet = model.grid[:,:,1].reshape((mr*mt))
		grid_points_r = grid_points_rm*np.sin(grid_points_tet)**2
		temp_points = model.temperature.reshape((mr*mt))
		plot.setYRange(0.9*np.amin(temp_points), 1.1*np.amax(temp_points))
		connect = np.full((mr*mt), 1); connect[mt-1::mt] = 0
		self.temp_line.setData(grid_points_r, temp_points, connect = connect)

class PopulationsWidget(QtGui.QWidget):
	def __init__(self, parent):
		self.parent = parent
		super(PopulationsWidget, self).__init__()
		layout = QtGui.QGridLayout()
		self.setLayout(layout)
		self.model_index = -1
		popul_plot = pg.GraphicsWindow()
		self.popul_plot = pg.PlotItem()
		self.popul_plot.setMenuEnabled(False)
		self.popul_plot.setMouseEnabled(False)
		# self.popul_plot.setAutoVisible(y = True)
		# self.popul_plot.enableAutoRange()
		
		self.popul_line = pg.PlotDataItem(pen = pg.mkPen(color = 'k', width = 2))
		self.pcolormesh_popul = pg.widgets.MatplotlibWidget.MatplotlibWidget()
		fig = self.pcolormesh_popul.getFigure()
		win = fig.canvas.window() 
		toolbar = win.findChild(QtGui.QToolBar)
		toolbar.setVisible(False)
		self.pcolormesh_popul_plot = fig.subplots()
		self.pcolormesh_popul_plot.axis('equal')
		self.pcolormesh_popul_plot.set_adjustable('box')
		self.popul_plot.addItem(self.popul_line)
		popul_plot.addItem(self.popul_plot)
		layout.addWidget(popul_plot, 0, 1, 5, 1)
		layout.addWidget(self.pcolormesh_popul, 0, 0, 5, 1)
		self.data_line = QtGui.QLabel('')
		self.prev_btn = QtGui.QPushButton('Previous')
		self.next_btn = QtGui.QPushButton('Next')
		self.log_btn = QtGui.QPushButton('Log scale')
		# name_label = QtGui.QLabel('File name')
		layout.addWidget(self.data_line, 0, 2)
		layout.addWidget(self.prev_btn, 1, 2)
		layout.addWidget(self.next_btn, 2, 2)
		layout.addWidget(self.log_btn, 3, 2)
		# layout.addWidget(name_label, 4, 1); layout.addWidget(self.name_text_line, 4, 2)
		layout.setColumnStretch(2, 0)
		layout.setColumnStretch(1, 2); layout.setColumnStretch(0, 3)
		# layout.setRowStretch(0, 0); layout.setRowStretch(1, 0); layout.setRowStretch(2, 0); layout.setRowStretch(3, 0)
		layout.setRowStretch(4, 10); layout.setRowStretch(0, 0)
		self.data_index = 0
		self.log_scale = False

	def load_model(self, model_index):
		try:
			self.next_btn.clicked.disconnect()
		except:
			pass
		try:
			self.prev_btn.clicked.disconnect()
		except:
			pass
		try:
			self.log_btn.clicked.disconnect()
		except:
			pass
		if model_index >= 0 and Main.models[model_index].populations_calculated:
			model = Main.models[model_index]
			self.model_index = model_index
			self.nh, self.ne, self.levels = model.get_populations()
			# if model.field = 'dipole':
			self.grid = model.grid
			self.data_index = 0
			self.headers = ['NH', 'Ne', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7']
			self.data_line.setText(self.headers[self.data_index])
			self.update_plots()
			self.next_btn.clicked.connect(self.increment_index)
			self.prev_btn.clicked.connect(self.decrement_index)
			self.log_btn.clicked.connect(self.set_scale)
		else:
			self.popul_line.setData([],[])
			self.data_index = -1
			self.update_plots()
			print('axes cleared')

	def set_scale(self):
		self.log_scale = not self.log_scale
		self.update_plots()

	def increment_index(self):
		self.data_index += 1
		self.data_index = self.data_index%9
		self.data_line.setText(self.headers[self.data_index])
		self.update_plots()

	def decrement_index(self):
		self.data_index -= 1
		self.data_index = self.data_index%9
		self.data_line.setText(self.headers[self.data_index])
		self.update_plots()

	def get_data(self):
		if self.data_index > 1:
			data = self.levels[:,:,self.data_index-2]
		else:
			data = [self.nh, self.ne][self.data_index]
		if self.log_scale:
			return np.ma.log10(data)
		else:
			return data

	def update_plots(self):
		def get_c(z):
			mx = len(z[:,0])
			my = len(z[0,:])
			c = np.zeros((mx,my))
			for i in range(mx-1):
				for j in range(my-1):
					c[i,j] = z[i,j] + z[i+1,j] + z[i,j+1] + z[i+1,j+1]       
			c = c/4
			return(c)
		# model = Main.models[self.model_index]
		if self.data_index >=0:
			plot = self.popul_plot
			mr, mt, dim = self.grid.shape
			grid_points_rm = self.grid[:,:,0].reshape((mr*mt))
			grid_points_tet = self.grid[:,:,1].reshape((mr*mt))
			grid_points_r = grid_points_rm*np.sin(grid_points_tet)**2
			grid_x = self.grid[:,:,0]*np.sin(self.grid[:,:,1])**3
			grid_y = self.grid[:,:,0]*np.sin(self.grid[:,:,1])**2*np.cos(self.grid[:,:,1])
			data = self.get_data()
			data_points = data.reshape((mr*mt))
			fig = self.pcolormesh_popul.getFigure()
			fig.clear()
			self.pcolormesh_popul_plot = fig.subplots()
			self.pcolormesh_popul_plot.axis('equal')
			self.pcolormesh_popul_plot.set_adjustable('box')
			data_map = self.pcolormesh_popul_plot.pcolormesh(grid_x,grid_y,get_c(data), cmap = 'Greys')
			cbar = fig.colorbar(data_map, ax = self.pcolormesh_popul_plot, format='%.1f')
			self.pcolormesh_popul.draw()
			plot.setYRange(0.9*np.amin(data_points), 1.1*np.amax(data_points))
			connect = np.full((mr*mt), 1); connect[mt-1::mt] = 0
			self.popul_line.setData(grid_points_r, data_points, connect = connect)
		else:
			self.popul_line.setData([], [])
			fig = self.pcolormesh_popul.getFigure()
			fig.clear()
			self.pcolormesh_popul_plot = fig.subplots()
			self.pcolormesh_popul_plot.axis('equal')
			self.pcolormesh_popul_plot.set_adjustable('box')
			data_map = self.pcolormesh_popul_plot.pcolormesh([[0,0],[1,1]],[[0,1],[0,1]],[[0]], cmap = 'Greys')
			cbar = fig.colorbar(data_map, ax = self.pcolormesh_popul_plot, format='%.1f')
	 
class ModelView(QtGui.QWidget):
	def __init__(self, parent):
		self.parent = parent
		super(ModelView, self).__init__()
		layout = QtGui.QGridLayout()
		self.setLayout(layout)
		self.tabs = QtGui.QTabWidget()
		layout.addWidget(self.tabs, 0, 0)
		self.model_grid_view = ModelGridView(self)
		self.model_temperature = TemperatureWidget(self)
		self.model_populations = PopulationsWidget(self)
		self.grid_view_tab = self.tabs.addTab(self.model_grid_view, 'Grid')
		self.temperature_tab = self.tabs.addTab(self.model_temperature, 'Temperature')
		self.populations_tab = self.tabs.addTab(self.model_populations, 'Populations')
		self.status_line = QtGui.QLabel('')
		layout.addWidget(self.status_line, 1, 0)
		self.model_index = -1
		self.tabs.tabBarClicked.connect(self.tab_clicked)

	def tab_clicked(self, index):
		# print(self.grid_updated)
		self.tabs.widget(index).load_model(self.model_index)


	def load_model(self, model_index):
		load_model_index = model_index
		if model_index >= 0:
			model = Main.models[model_index]
			if model.field.type == 'dipole':
				self.set_status_line('Loaded model '+model.name, content = 'msg')
			else:
				self.set_status_line('Field type '+model.field.type+" isn't supported yet", content = 'err')
				load_model_index = -1
		self.model_index = load_model_index
		self.model_grid_view.load_model(self.model_index)
		self.tabs.setCurrentIndex(self.grid_view_tab)
		
		# self.model_temperature.load_model(load_model_index)

	def set_status_line(self, message, content = 'msg'):
		contents_style = {'msg' : "color : #000000",
											'err' : "color : #aa0000"}
		self.status_line.setText(message)
		self.status_line.setStyleSheet(contents_style[content])

	def reset(self):
		self.load_model(-1)


class MainWindow(QtGui.QWidget):
	def __init__(self):
		super(MainWindow, self).__init__()
		layout = QtGui.QVBoxLayout()
		self.setLayout(layout)
		self.model_manager = ModelManager(self)
		layout.addWidget(self.model_manager)
		self.model_view = ModelView(self)
		layout.addWidget(self.model_view)
		# layout.setStretch(0,0)
		layout.setStretch(1,1)



class Main:
	models = []

	def __init__(self):
		self.app = QtGui.QApplication([])
		self.widget = MainWindow()
		# layout = QtGui.QGridLayout()
		# self.widget.setLayout(layout)
		# self.widget.add_model_window = QtGui.QWidget()
		# self.widget.add_model_widget.setLayout(QtGui.QGridLayout())
		# create_model_btn = QtGui.QPushButton('Add model')
		# self.widget.create_model_btn.clicked.connect(self.add_model_window)
		# layout.addWidget(create_model_btn, 0, 0)

	def add_model(self, model):
		self.models.append(model)
		self.widget.model_list.addItem(model.name)

	def main(self):
		def create_directories(local_path):
			created_dir = os.getcwd() + '/'
			for directory in local_path.split('/')[:-1]:
				created_dir = created_dir + directory + '/' 
				try:
					os.mkdir(created_dir)
				except FileExistsError:
					pass
		create_directories(model_data_dir)
		create_directories(populations_dir)
		create_directories(hotspot_data_dir)
		self.widget.show()
		self.app.exec_()



main = Main()
main.main()

# app = QtGui.QApplication([])

# model = Model('hart94_7_75_new', new = True)
# # print(model.field.type)
# ## Always start by initializing Qt (only once per application)

# ## Define a top-level widget to hold everything
# w = QtGui.QWidget()

# create_model_btn = QtGui.QPushButton('Create model')
# create_model_btn.clicked.connect()
# layout.addWidget(create_model_btn, 0, 0)   # button goes in upper-left

# btn = QtGui.QPushButton('clear')
# listw = QtGui.QListWidget()
# plot = pg.widgets.MatplotlibWidget.MatplotlibWidget()# PlotWidget()

# def clear_list():
#   listw.clear()

# def click(event):
#   xd = event.xdata
#   yd = event.ydata
#   borders = prof.init_one_border(yd, xd)
#   try:
#     listw.addItem(('{:.3f}, '*8 + '{:.3f}').format(*borders))
#   except:
#     pass

# ## Create some widgets to be placed inside

# btn.clicked.connect(clear_list)
# ## Create a grid layout to manage the widgets size and position
# layout = QtGui.QGridLayout()
# w.setLayout(layout)

# ## Add widgets to the layout in their proper positions
# layout.addWidget(btn, 0, 0)   # button goes in upper-left
# layout.addWidget(listw, 1, 0)  # list widget goes in bottom-left
# layout.addWidget(plot, 0, 1, 3, 1)  # plot goes on right side, spanning 3 rows
# layout.setColumnStretch(1,100)

# m=50

# i=90; psi=90; alpha=15

# prof.init_orientation(i, psi, alpha)
# prof.init_field('dipole', 1e-7, 2.2, 3.0, 1.0, 1.2)
# field_borders, grid, dS = prof.init_field_borders(m)

# x = grid[1,:,:]
# y = grid[0,:,:]
# z = field_borders[0,:,:]

# def get_c(z):
#   mx = len(z[:,0])
#   my = len(z[0,:])
#   c = np.zeros((mx,my))
#   for i in range(mx-1):
#     for j in range(my-1):
#       c[i,j] = z[i,j] + z[i+1,j] + z[i,j+1] + z[i+1,j+1]       
#   c = c/4
#   return(c)

# c = get_c(z)

# fig = plot.getFigure()
# win = fig.canvas.window()
# toolbar = win.findChild(QtGui.QToolBar)
# toolbar.setVisible(False)

# subplot = fig.subplots()
# subplot.axis('equal')
# subplot.set_adjustable('box')
# testmap = subplot.pcolormesh(x,y,c, cmap = 'Greys')
# cbar = fig.colorbar(testmap, ax = subplot, format='%.1f')
# fig.canvas.mpl_connect("button_press_event", click)
# plot.draw()

# ## Display the widget as a new window
# w.show()

# ## Start the Qt event loop
# app.exec_()
