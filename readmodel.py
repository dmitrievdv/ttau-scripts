import numpy as np

global_parameters = {"Mdot" : -1.0,
          "Mstar" : -1.0, 
          "Tstar" : -1.0, 
          "Rstar" : -1.0,
          "Tmax" : 8000,
          "field_type" : "none",
          "first_border" : -1.0,
          "second_border" : -1.0,
          "in_cut" : -1.0,
          "out_cut" : -1.0,
          "hot_spot" : False,
          "Thot" : -1.0,
          "Dhot" : -1.0,
          "populations" : ''}

global_parameters_new = {"Mdot" : -1.0,
          "Mstar" : -1.0, 
          "Tstar" : -1.0, 
          "Rstar" : -1.0,
          "Tmax" : 8000,
          "field_type" : "none",
          "first_border" : -1.0,
          "second_border" : -1.0,
          "in_cut" : -1.0,
          "out_cut" : -1.0,
          "hotspot" : '',
          "Thot" : -1.0,
          "Dhot" : -1.0,
          "populations" : ''}

parameter_names = {"Mdot" : "Mdot",
             "Mstar" : "Mstar", 
             "Tstar" : "Tstar", 
             "Rstar" : "Rstar",
             "Tmag" : "Tmax", 
             "field type" : "field_type",
             "first" : "first_border",
             "second" : "second_border",
             "inner" : "in_cut",
             "outer" : "out_cut",
             "hot spot:" : "hotspot",
             "hot spot" : "hotspot",
             "Thot" : "Thot",
             "Dhot" : "Dhot",
             "populations" : "populations"}

parameter_read_funcs = {"Mdot" : float,
          "Mstar" : float, 
          "Tstar" : float, 
          "Rstar" : float,
          "Tmag" : float,
          "field type" : lambda x: x,
          "first" : float,
          "second" : float,
          "inner" : float,
          "outer" : float,
          "hot spot:" : lambda x: True,
          "hot spot" : str,
          "Thot" : float,
          "Dhot" : float,
          "populations" : str}

important_parameters = {"Mdot" : False,
          "Mstar" : False, 
          "Tstar" : False, 
          "Rstar" : False,
          "Tmax" : False,
          "field_type" : False,
          "first_border" : False,
          "second_border" : False}


class Error(Exception):
  pass

class NoSuchModelError(Error):
  def __init__(self, model_name):
    self.model_name = model_name

class NoDataForLevel(Error):
  def __init__(self, model_name, level):
    self.model_name = model_name
    self.level = level

class ImportantParameterNotSet(Error):
  def __init__(self, parameter, parameter_key):
    self.parameter = parameter
    self.parameter_key = parameter_key

class InvalidParameterValue(Error):
  def __init__(self, parameter_key):
    self.parameter_key = parameter_key

class UpLevelLowerThenLowLevel(Error):
  pass

def read_parameters_from_string(model_name, data):
  parameters = global_parameters_new
  parameters['populations'] = model_name
  parameters['hotspot'] = ''
  print(global_parameters_new['hotspot'])
  important_parameters_set_status = {"Mdot" : False,
          "Mstar" : False, 
          "Tstar" : False, 
          "Rstar" : False,
          "Tmax" : False,
          "field_type" : False,
          "first_border" : False,
          "second_border" : False}
  # print(important_parameters_set_status)
  data_lines = data.split('\n')
  for line in data_lines:
    if line != '':
      content, comment = tuple(line.split('#'))
      comment = comment.strip()
      content = content.strip()
      try:
        parameters[parameter_names[comment]] = parameter_read_funcs[comment](content)
      except KeyError:
        pass
      except ValueError:
        raise InvalidParameterValue(parameter_names[comment])
      try:
        important_parameters_set_status[parameter_names[comment]] = True
      except KeyError:
        pass

  if(parameters['out_cut'] == -1.0):
    parameters['out_cut'] = parameters['second_border']
  if(parameters['in_cut'] == -1.0):
    parameters['in_cut'] = 1.0

  parameter_keys = {v : k for k,v in parameter_names.items()}

  if(parameters['field_type'] == 'cone'):
    important_parameters_set_status['Tmax'] = True

  for key, value in important_parameters_set_status.items():
    if not value:
      raise ImportantParameterNotSet(key, parameter_keys[key])


  return parameters

def read_parameters(model_name, model_dir = 'models/data'):

  parameters = global_parameters
  parameters['populations'] = model_name
  parameters_file_name = model_name + '_data.dat'
  try:
    print(parameters_file_name)
    parameters_file = open(model_dir + '/'+parameters_file_name, 'r')
  except FileNotFoundError:
    raise NoSuchModelError(model_name)

  for line in parameters_file.readlines():
    # print(line)
    content, comment = tuple(line.split('#'))
    comment = comment.strip()
    content = content.strip()
    # print(comment, content)
    try:
      parameters[parameter_names[comment]] = parameter_read_funcs[comment](content)
    except KeyError:
      pass

  if(parameters['out_cut'] == -1.0): 
    parameters['out_cut'] = parameters['second_border']
  if(parameters['in_cut'] == -1.0): 
    parameters['in_cut'] = 1.0
  return parameters

def read_populations_file(model_name, u, l, field = 'dipole', model_dir = 'models/popul'):
  if(u < l):
    raise UpLevelLowerThenLowLevel
  try:
    populations_file = open(model_dir+'/'+model_name+'_popul.dat', 'r')
  except FileNotFoundError:
    raise NoSuchModelError(model_name)
  lines = populations_file.readlines()
  grid_precision_data = lines[0].split('#')[1].strip()
  if(field == 'dipole') :
    rm_precision = int(grid_precision_data.split()[0])
    dtet_precision = int(grid_precision_data.split()[1])
    try:
      maxlevel = int(lines[1].split()[5][-1])
    except:
      maxlevel = 2
    if(u > maxlevel):
      raise NoDataForLevel(model_name, u)
    get_column = {'rm' : lambda: 0, 'tet' : lambda: 1, 'Te' : lambda: 2,
            'nh' : lambda: 3, 'ne' : lambda: -1, 'n' : lambda x: 3+x}
    Te = np.zeros((rm_precision, dtet_precision))
    n_u = np.zeros((rm_precision, dtet_precision))
    n_l = np.zeros((rm_precision, dtet_precision))
    nh = np.zeros((rm_precision, dtet_precision))
    ne = np.zeros((rm_precision, dtet_precision))  
    grid = np.zeros((rm_precision, dtet_precision, 2))
    rm_position = 0
    tet_position = 0
    for line in lines[2:]:
      if(line.strip() == ''):
        tet_position = 0
        rm_position += 1
      else:
        columns = line.split()
        grid[rm_position,tet_position,0] = float(columns[get_column['rm']()])
        grid[rm_position,tet_position,1] = float(columns[get_column['tet']()])
        Te[rm_position,tet_position] = float(columns[get_column['Te']()])
        ne[rm_position,tet_position] = float(columns[get_column['ne']()])
        nh[rm_position,tet_position] = float(columns[get_column['nh']()])
        n_u[rm_position,tet_position] = float(columns[get_column['n'](u)])
        n_l[rm_position,tet_position] = float(columns[get_column['n'](l)])
        tet_position += 1
  elif(field == 'cone'):
    r_precision = int(grid_precision_data.split()[0])
    Te = np.zeros((r_precision))
    n_u = np.zeros((r_precision))
    n_l = np.zeros((r_precision))
    nh = np.zeros((r_precision))
    ne = np.zeros((r_precision))  
    grid = np.zeros((r_precision))
    get_column = {'r' : lambda: 0, 'Te' : lambda: 1,
            'nh' : lambda: 3, 'ne' : lambda: -1, 'n' : lambda x: 1+x}
    r_position = 0
    for line in lines[2:]:
      columns = line.split()
      grid[r_position] = float(columns[get_column['r']()])
      Te[r_position] = float(columns[get_column['Te']()])
      ne[r_position] = float(columns[get_column['ne']()])
      nh[r_position] = ne[r_position] + float(columns[get_column['nh']()])
      n_u[r_position] = float(columns[get_column['n'](u)])
      n_l[r_position] = float(columns[get_column['n'](l)])
      r_position += 1
    grid = grid/grid[0]


  return grid, Te, nh, ne, n_u, n_l 


def gen_data_file(model_name, file_content, model_dir = 'models'):
  data_file = open(model_dir + '/' + model_name, 'w')
  for comment, content in zip(file_content.keys(), file_content.values()):
    data_file.write(str(content) + ' # ' + str(comment) + '\n')