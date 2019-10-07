# These functions are intended to handle and analyze de novo NTF2 blueprints.
# Blueprint library is from Javier Castellanos.
# Benjamin Basanta - June 5th 2017 - Baker Lab - Bichemistry dept. - UW.

from Blueprint import *
import random
import math
from pyrosetta import *
from rosetta import *
import numpy as np
from glob import glob
import pandas as pnd
import json

def get_base_width(bp):
	# Base widths:
	# 3: 4x3
	# 5: 4x5
	# 7: 4x7
	segments = [i for i in bp.segment_dict.keys()]
	strands = [ i for i in segments if i[0]=='E']
	B1_A_pos = [ i[0] for i in bp.segment_dict['E3'].bp_data if i[2]=='EA'][0]
	E3_len = len(bp.segment_dict['E3'].bp_data)
	E4_len = len(bp.segment_dict['E4'].bp_data)

	last_pos_E3 = bp.segment_dict['E3'].bp_data[-1][0]
	B2_A_pos = [ i[0] for i in bp.segment_dict['E6'].bp_data if i[2]=='EA'][0]

	E6_init = bp.segment_dict['E6'].bp_data[0][0]
	E4_E5_RS = int([i[3] for i in bp.sspairs if i[0]=='4' and i[1]=='5'][0])
	E5_len = len(bp.segment_dict['E5'].bp_data)
	#E4_E5_delta = ( E4_E5_RS + E5_len ) - E4_len
	E6_N_to_Bul = B2_A_pos - E6_init

	inter_bulge_dist = (last_pos_E3-B1_A_pos) + (E4_E5_RS + E6_N_to_Bul )
	return inter_bulge_dist

def get_secstruct_from_PDB_info(fname):
	handle = open(fname,'r')
	lines = handle.readlines()
	handle.close()
	for line in lines:
		if len(line.split()) == 2:
			if line.split()[0] == 'SECSTRUCT':
				return line.split()[1]

def get_dict_line_from_PDB_info(fname,tag):
	handle = open(fname,'r')
	lines = handle.readlines()
	handle.close()
	for line in lines:
		if len(line.split()) > 2:
			if line.split()[0] == tag:
				return (' ').join(line.split()[1:])

def is_4x5(bp):
	if get_base_width(bp) == 5: return True
	else: return False

def CircularHBondConstraints(donor,acceptor): # return string for cst file
	hb_ang = np.deg2rad(180.)
	#hb_ang_tol=np.deg2rad(20.0) -> original
	hb_ang_tol=np.deg2rad(90.0)
	st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(donor,acceptor)
	#st = "AtomPair N %i O %i BOUNDED 2.7 3.3 0.2 \n" %(donor,acceptor)
	st+= "Angle N %i H %i O %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
	st+= "Angle H %i O %i C %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,acceptor,acceptor,hb_ang,hb_ang_tol)
	return st

def PerfectHelixCst(init_blue,helixn):
	blue = Blueprint(data=[ [ x for x in i ] for i in init_blue.bp_data ] )
	blue.reindex_blueprint(start=1)
	hnum = 'H'+str(helixn)
	H = blue.segment_dict[hnum]
	HN = int(H.bp_data[0][0])
	HC = int(H.bp_data[-1][0])
	cst = ''
	phi = np.deg2rad(-63)
	psi = np.deg2rad(-42)
	sd = np.deg2rad(50)
	cst += "Dihedral N %i CA %i C %i N %i HARMONIC %0.2f %0.2f \n" %(HN,HN,HN,HN+1,psi,sd) #psi
	for pos in range(HN+1,HC-1):
		cst += "Dihedral C %i N %i CA %i C %i HARMONIC %0.2f %0.2f \n" %(pos-1,pos,pos,pos,phi,sd) #phi
		cst += "Dihedral N %i CA %i C %i N %i HARMONIC %0.2f %0.2f \n" %(pos,pos,pos,pos+1,psi,sd) #psi
	cst += "Dihedral C %i N %i CA %i C %i HARMONIC %0.2f %0.2f \n" %(HC-1,HC,HC,HC,phi,sd) #phi
	return cst

class NTF2_sheet():

	sheet_atributes_Dict = {
	"base_width":[3,5],
	"long_arm_l":[2,3,4],
	"short_arm_l":[1,2],
	"Second_bulge_E3":[True,False],
	"Second_b_place":[1,2],
	"E3_MainBulgeCurve":[1,2,3],
	"E3_SecBulgeCurve":[1,2,3],
	"ExtendedE6":[True,False],
	"ExtendedE4":[True,False],
	"CurvedLongArm":[True,False]
	}

	def validate_sheet(self):
		exit_status_OK = True
		if (type(self.base_width) is not int) or \
		(type(self.long_arm_l) is not int) or \
		(type(self.short_arm_l) is not int) or \
		(type(self.Second_bulge_E3) is not bool) or \
		(type(self.E3_MainBulgeCurve) is not int) or\
		(type(self.CurvedLongArm) is not bool ):
			raise ValueError("First check not passed, values are the wrong type")
			exit_status_OK = False
		if self.Second_bulge_E3:
			if (type(self.Second_b_place) is not int) or \
			(type(self.E3_SecBulgeCurve) is not int):
				raise Exception("Incompatible values")
		if (self.E3_MainBulgeCurve == 1) and (type(self.ExtendedE6) is not bool):
			raise Exception("Incompatible values")
		if not exit_status_OK: raise ValueError("You asked for a second bulge in E3 but it's placement is not an integer")
		# Now check to make sure that the selected values are allowed:
		if (self.base_width not in NTF2_sheet.sheet_atributes_Dict["base_width"]):
			raise Exception("Incompatible values")
		if (self.long_arm_l not in NTF2_sheet.sheet_atributes_Dict["long_arm_l"]):
			raise Exception("Incompatible values")
		if (self.short_arm_l not in NTF2_sheet.sheet_atributes_Dict["short_arm_l"]):
			raise Exception("Incompatible values")
		if (self.Second_bulge_E3 not in NTF2_sheet.sheet_atributes_Dict["Second_bulge_E3"]):
			raise Exception("Incompatible values")
		if (self.E3_MainBulgeCurve not in NTF2_sheet.sheet_atributes_Dict["E3_MainBulgeCurve"]):
			#exit_status_OK = False
			raise Exception("Incompatible values")
		if self.Second_bulge_E3:
			if (self.Second_b_place not in NTF2_sheet.sheet_atributes_Dict["Second_b_place"]) or \
			(self.E3_SecBulgeCurve not in NTF2_sheet.sheet_atributes_Dict["E3_SecBulgeCurve"]):
				exit_status_OK = False
		if (self.E3_MainBulgeCurve == 1) and (self.ExtendedE6 not in NTF2_sheet.sheet_atributes_Dict["ExtendedE6"]):
			exit_status_OK = False
		if not exit_status_OK: raise ValueError("One of the input values is not allowed, check sheet_atributes_dict")
		# Now some basic checks:
		if (self.long_arm_l == 2) and self.Second_bulge_E3:
			raise ValueError("Can't have a second bulge on E3 if long_arm_l=2")
			exit_status_OK = False
		if self.Second_bulge_E3 and (self.Second_b_place == 2) and (self.long_arm_l <= 3):
			raise ValueError("Second bulge pos is 2 but long_arm_l<=3")
			exit_status_OK = False
		if (self.base_width > 3) and (self.long_arm_l == 4) and not self.Second_bulge_E3:
			raise ValueError("Sheet is too long, make it shorter or add 2nd bulge on E3")
			exit_status_OK = False
		if (self.E3_MainBulgeCurve > 1) and self.ExtendedE6:
			raise ValueError("The sheet is too curved to allow extention of E6")
			exit_status_OK = False
		# Extended E4 and base_width = 3 checks:
		if not self.ExtendedE4:
			if (self.base_width == 3):
				raise ValueError("Base width 3 needs, short_arm_l == 2 and extended E4")
				exit_status_OK = False
		if (self.short_arm_l == 1):
			if self.ExtendedE4:
				raise ValueError("short_arm_l == 1 is incompatible with E4 extention")
				exit_status_OK = False
		if (self.short_arm_l == 1):
			if (self.base_width == 3):
				raise ValueError("Base width 3 needs, short_arm_l == 2 and extended E4")
				exit_status_OK = False
		if (self.short_arm_l == 2) and not self.ExtendedE4:
			raise ValueError("short_arm_l = 2 requires extended E4")
			exit_status_OK = False
		# Curvature and long_arm_length checks:
		if (self.long_arm_l > 2) and not self.Second_bulge_E3 and (self.E3_MainBulgeCurve < 2):
			raise ValueError("For long_arm_l > 2, if there is no secondary bulge curve at main bulge must be higher")
			exit_status_OK = False
		if (self.long_arm_l > 3) and not self.Second_bulge_E3 and (self.E3_MainBulgeCurve < 3):
			raise ValueError("For long_arm_l > 3, if there is no secondary bulge curve at main bulge must be 3")
			exit_status_OK = False
		if (self.long_arm_l >= 3) and (self.E3_MainBulgeCurve < 2):
			raise ValueError("Sheet is too extended")
			exit_status_OK = False
		if (self.long_arm_l > 3) and self.Second_bulge_E3 and (self.Second_b_place == 1) and (self.E3_SecBulgeCurve != 2):
			raise ValueError("When long arm length = 4, there is a second bulge and that bulge is close to the main bulge, then the curvature at sec bulge must be = 2")
			exit_status_OK = False
		if self.CurvedLongArm and self.Second_bulge_E3:
			raise ValueError("Incompatible CurvedLongArm and Second E3 bulge")
			exit_status_OK = False
		if self.CurvedLongArm and (self.long_arm_l < 3):
			raise ValueError("The long arm is too short to curve it.")
			exit_status_OK = False
		if exit_status_OK:
			self.sheet_data["Validated"] = True
		else:
			self.sheet_data["Validated"] = False

	def soft_reset_sheet(self):
		self.base_width = self.saved_input["base_width"]
		self.long_arm_l = self.saved_input["long_arm_l"]
		self.short_arm_l = self.saved_input["short_arm_l"] # Chort arm length, 2*X residues.
		self.Second_bulge_E3 = self.saved_input["Second_bulge_E3"] # Is there a second bulge at E3? boolean
		self.Second_b_place = self.saved_input["Second_b_place"] # If there is a second bulge at E3, where in E3? places are 2, 4 or 6 places before the main E3 bulge.
		self.E3_MainBulgeCurve = self.saved_input["E3_MainBulgeCurve"] # What is the overall curvature of E3, E4 and E5 after E3 bulge? If there is a second bulge in E3 this parameter is not used. 1: No curvature, 2 and 3 increasing degrees of curvature.
		self.E3_SecBulgeCurve = self.saved_input["E3_SecBulgeCurve"]
		self.ExtendedE6 = self.saved_input["ExtendedE6"]
		self.ExtendedE4 = self.saved_input["ExtendedE4"]
		self.CurvedLongArm = self.saved_input["CurvedLongArm"]

		self.csts = None
		self.pairings = None
		self.blueprint = None

	def update_sheet_data(self):
		self.sheet_data["base_width"] = self.base_width
		self.sheet_data["long_arm_l"] = self.long_arm_l
		self.sheet_data["short_arm_l"] = self.short_arm_l
		self.sheet_data["Second_bulge_E3"] = self.Second_bulge_E3
		self.sheet_data["Second_b_place"] = self.Second_b_place
		self.sheet_data["E3_MainBulgeCurve"] = self.E3_MainBulgeCurve
		self.sheet_data["E3_SecBulgeCurve"] = self.E3_SecBulgeCurve
		self.sheet_data["ExtendedE6"] = self.ExtendedE6
		self.sheet_data["ExtendedE4"] = self.ExtendedE4
		self.sheet_data["CurvedLongArm"] = self.CurvedLongArm

	def create_sheet(self,validate=True):
		"""
		Creates a sheet with input parameters, parameters not specidied are randomized
		This function also validates the combination of parameters and creates blueprint and constraint lines.
		"""
		if not self.base_width:
			print( "Base width not set, randomizing")
			self.base_width = random.choice(NTF2_sheet.sheet_atributes_Dict["base_width"])
		if (self.base_width == 3) and (self.short_arm_l is None) and (self.ExtendedE4 is None):
			print( "Base width 3 needs short_arm_l=2 and estended E4, the latter are not defined, so I'm setting them to the required values")
			self.short_arm_l = 2
			self.ExtendedE4 = True
		if not self.long_arm_l:
			print( "Long arm length not set, randomizing")
			self.long_arm_l = random.choice(NTF2_sheet.sheet_atributes_Dict["long_arm_l"])
		if not self.short_arm_l:
			print( "Short arm length not set, randomizing")
			self.short_arm_l = random.choice(NTF2_sheet.sheet_atributes_Dict["short_arm_l"])
		if (self.Second_bulge_E3 is None):
			print( "Existence of secondary bulge in E3 not set, randomizing")
			self.Second_bulge_E3 = random.choice(NTF2_sheet.sheet_atributes_Dict["Second_bulge_E3"])
		if not self.E3_MainBulgeCurve:
			print( "Curvature after E3 main bulge not set, randomizing")
			self.E3_MainBulgeCurve = random.choice(NTF2_sheet.sheet_atributes_Dict["E3_MainBulgeCurve"])
		if self.Second_bulge_E3 and not self.E3_SecBulgeCurve:
			print( "Curvature after E3 secondary bulge not set, randomizing")
			self.E3_SecBulgeCurve = random.choice(NTF2_sheet.sheet_atributes_Dict["E3_SecBulgeCurve"])
		if self.Second_bulge_E3 and not self.Second_b_place :
			print( "Placement of secondary bulge in E3 not set, randomizing")
			self.Second_b_place = random.choice(NTF2_sheet.sheet_atributes_Dict["Second_b_place"])
		if (self.ExtendedE6 is None) and (self.E3_MainBulgeCurve == 1):
			print( "Extension of E6 not set and curvature on E3 low enough to allow it, randomizing")
			self.ExtendedE6 = random.choice(NTF2_sheet.sheet_atributes_Dict["ExtendedE6"])
		if (self.ExtendedE4 is None) and (self.short_arm_l == 2):
			print( "Extended E4 will be set to true because self.short_arm_l = 2")
			self.ExtendedE4 = True
		if (self.CurvedLongArm is None) and not self.Second_bulge_E3 and (self.long_arm_l > 2):
			print( "Long arm length allows for curvature without bulge, randomizing")
			self.CurvedLongArm = random.choice(NTF2_sheet.sheet_atributes_Dict["CurvedLongArm"])

		if validate:
			try:
				self.validate_sheet()
			except ValueError as e:
				print(e)
				print( "Initialization failed, some parameters are incompatible, trying again, but check the problem is not in the input values")
				self.soft_reset_sheet()
				self.create_sheet(validate=True)

			self.update_sheet_data()
			self.blueprint = self.produce_bp()
			self.csts = self.produce_csts()
		#else:
			#print("Sheet won't be validated, remember to update_sheet_data(), self.blueprint = self.produce_bp() and self.csts = self.produce_csts() once sheet is validated")

	def __init__(self, arch_dist=None, base_width = None, long_arm_l = None, short_arm_l = None, Second_bulge_E3 = None,\
	Second_b_place = None, E3_MainBulgeCurve = None, E3_SecBulgeCurve = None, ExtendedE6 = None, ExtendedE4 = None, validate_sheet=True, CurvedLongArm=False):
		#print("Creating NTF2")
		self.saved_input={
		"base_width" : base_width, \
		"long_arm_l" : long_arm_l, \
		"short_arm_l" : short_arm_l, \
		"Second_bulge_E3" : Second_bulge_E3, \
		"Second_b_place" : Second_b_place, \
		"E3_MainBulgeCurve" : E3_MainBulgeCurve, \
		"E3_SecBulgeCurve" : E3_SecBulgeCurve, \
		"ExtendedE6" : ExtendedE6, \
		"ExtendedE4" : ExtendedE4, \
		"CurvedLongArm" : CurvedLongArm \
		}

		self.sheet_data = {}
		self.base_width = base_width # Relative separation between E3 bulge and E5 bulge
		self.long_arm_l = long_arm_l # Length of long arm: 1+2*X aas, or 2*(X+1) if there is an additional bulge
		self.short_arm_l = short_arm_l # Chort arm length, 2*X residues.
		self.Second_bulge_E3 = Second_bulge_E3 # Is there a second bulge at E3? boolean
		self.Second_b_place = Second_b_place # If there is a second bulge at E3, where in E3? places are 2, 4 or 6 places before the main E3 bulge.
		self.E3_MainBulgeCurve = E3_MainBulgeCurve # What is the overall curvature of E3, E4 and E5 after E3 bulge? If there is a second bulge in E3 this parameter is not used. 1: No curvature, 2 and 3 increasing degrees of curvature.
		self.E3_SecBulgeCurve = E3_SecBulgeCurve
		self.ExtendedE6 = ExtendedE6
		self.ExtendedE4 = ExtendedE4
		self.CurvedLongArm = CurvedLongArm

		self.pairings = None
		self.blueprint = None
		self.csts = None

		self.arch_dist = arch_dist
		if validate_sheet:
			self.create_sheet(validate=True)
		else:
			self.create_sheet(validate=False)

	def create_basic_sheet(self): # This is a quick way to create a simple, validated sheet structure.
		self.base_width = 5
		self.long_arm_l = 2
		self.short_arm_l = 1
		self.Second_bulge_E3 = False
		self.Second_b_place = None
		self.E3_MainBulgeCurve = 1
		self.E3_SecBulgeCurve = None
		self.ExtendedE6 = False
		self.ExtendedE4 = False
		self.CurvedLongArm = False

		self.pairings = None
		self.blueprint = None
		self.csts = None
		self.create_sheet(validate=True)

	def create_random_sheet(self):
		self.base_width = None
		self.long_arm_l = None
		self.short_arm_l = None
		self.Second_bulge_E3 = None
		self.Second_b_place = None
		self.E3_MainBulgeCurve = None
		self.E3_SecBulgeCurve = None
		self.ExtendedE6 = None
		self.ExtendedE4 = None
		self.CurvedLongArm = None

		self.pairings = None
		self.blueprint = None
		self.csts = None
		self.create_sheet()

	def hard_reset_sheet(self):
		self.saved_input={
		"base_width" : None, \
		"long_arm_l" : None, \
		"short_arm_l" : None, \
		"Second_bulge_E3" : None, \
		"Second_b_place" : None, \
		"E3_MainBulgeCurve" : None, \
		"E3_SecBulgeCurve" : None, \
		"ExtendedE6" : None, \
		"ExtendedE4" : None, \
		"CurvedLongArm" : None \
		}

		self.base_width = None
		self.long_arm_l = None
		self.short_arm_l = None
		self.Second_bulge_E3 = None
		self.Second_b_place = None
		self.E3_MainBulgeCurve = None
		self.E3_SecBulgeCurve = None
		self.ExtendedE6 = None
		self.ExtendedE4 = None
		self.CurvedLongArm = None

		self.csts = None
		self.pairings = None
		self.blueprint = None

	def produce_csts(self):
		# Arch distance constaint limit: 27A between E6 bulge and E3 N:
		E3N = self.blueprint.segment_dict['E1'].bp_data[0][0]
		E6_bulge_INDEX_from_N = 2*(self.short_arm_l)
		E6 = self.blueprint.segment_dict['E4']
		E6Bulge_A = E6.bp_data[E6_bulge_INDEX_from_N][0]
		CST_Longest_arch_dist = "AtomPair CA %i CA %i FLAT_HARMONIC 20.0 3.0 10.0"%(E3N,E6Bulge_A)
		if self.arch_dist:
			CST_Longest_arch_dist = "AtomPair CA %i CA %i HARMONIC %0.2f 3.0"%(E3N,E6Bulge_A,self.arch_dist)
		# E3 main curvature
		if self.ExtendedE4:
			relative_MainE4_curve_center = self.base_width + 2
		else:
			relative_MainE4_curve_center = self.base_width
		MainE4_curve_center = self.blueprint.segment_dict['E2'].bp_data[relative_MainE4_curve_center-1][0]
		MainE4_curve_root_1 = self.blueprint.segment_dict['E2'].bp_data[0][0]
		MainE4_curve_root_2 = MainE4_curve_center + 2
		#"E3_MainBulgeCurve":[1,2,3],
		mean_angle_value = math.radians(165 - self.E3_MainBulgeCurve*10)
		std_dev = math.radians(5)
		#
		# Add lines for actual cst string
		CST_line_MainE4 = "Angle CA %d CA %d CA %d HARMONIC %0.2f %0.2f" \
		%(MainE4_curve_root_1,MainE4_curve_center,MainE4_curve_root_2,mean_angle_value,std_dev)


		E4_len = len(self.blueprint.segment_dict['E2'].bp_data)
		relative_MainE5_curve_center = E4_len - relative_MainE4_curve_center - 2
		MainE5_curve_center = self.blueprint.segment_dict['E3'].bp_data[relative_MainE5_curve_center][0]
		MainE5_curve_root_1 = MainE5_curve_center - 2
		MainE5_curve_root_2 = MainE5_curve_center + 2

		#
		# Add lines for actual cst string
		CST_line_MainE5 = "Angle CA %d CA %d CA %d HARMONIC %0.2f %0.2f" \
		%(MainE5_curve_root_1,MainE5_curve_center,MainE5_curve_root_2,mean_angle_value,std_dev)
		#

		CST_RETURN_LIST = [CST_line_MainE4, CST_line_MainE5, CST_Longest_arch_dist]

		# E3 secondary curvature if applicable:
		if self.Second_bulge_E3:
			mean_angle_value = math.radians(165 - self.E3_SecBulgeCurve*10)
			std_dev = math.radians(5)
			# "E3_SecBulgeCurve":[1,2,3],
			relative_SecE4_curve_center = relative_MainE4_curve_center + 2*self.Second_b_place + 2
			SecE4_curve_center = self.blueprint.segment_dict['E2'].bp_data[relative_SecE4_curve_center-1][0]
			SecE4_curve_root_1 = SecE4_curve_center - 2
			SecE4_curve_root_2 = SecE4_curve_center + 2

			#
			# Add lines for actual cst string
			CST_line_SecE4 = "Angle CA %d CA %d CA %d HARMONIC %0.2f %0.2f" \
			%(SecE4_curve_root_1,SecE4_curve_center,SecE4_curve_root_2,mean_angle_value,std_dev)
			#
			CST_RETURN_LIST.append(CST_line_SecE4)

			#if (E4_len - relative_SecE4_curve_center) > 3:
			relative_SecE5_curve_center = E4_len - relative_SecE4_curve_center # keep at the same level?? Yes, works better.
			SecE5_curve_center = self.blueprint.segment_dict['E3'].bp_data[relative_SecE5_curve_center][0]
			SecE5_curve_root_1 = SecE5_curve_center - 2
			SecE5_curve_root_2 = SecE5_curve_center + 2
			CST_line_SecE5 = "Angle CA %d CA %d CA %d HARMONIC %0.2f %0.2f" \
			%(SecE5_curve_root_1,SecE5_curve_center,SecE5_curve_root_2,mean_angle_value,std_dev)

			CST_RETURN_LIST.append(CST_line_SecE5)
		###################################

		elif self.CurvedLongArm: # The position of the secondary curvature is fixed!
			mean_angle_value = math.radians(150) #155
			std_dev = math.radians(5)
			# Curve on E3:
			relative_SecE3_curve_center = -1*(relative_MainE4_curve_center + 1) - 5 # -4
			SecE3_curve_center = self.blueprint.segment_dict['E1'].bp_data[relative_SecE3_curve_center][0]
			#SecE3_curve_root_1 = SecE3_curve_center - 2
			SecE3_curve_root_1 = self.blueprint.segment_dict['E1'].bp_data[0][0]
			#SecE3_curve_root_2 = SecE3_curve_center + 2
			SecE3_curve_root_2 = SecE3_curve_center + 4


			CST_line_SecE3 = "Angle CA %d CA %d CA %d HARMONIC %0.2f %0.2f" \
			%(SecE3_curve_root_1,SecE3_curve_center,SecE3_curve_root_2,mean_angle_value,std_dev)
			#
			CST_RETURN_LIST.append(CST_line_SecE3)

			# Now for E4:
			relative_SecE4_curve_center = relative_MainE4_curve_center + 5 # +4
			SecE4_curve_center = self.blueprint.segment_dict['E2'].bp_data[relative_SecE4_curve_center-1][0]
			#
			#SecE4_curve_root_1 = SecE4_curve_center - 2
			SecE4_curve_root_1 = SecE4_curve_center - 4
			#SecE4_curve_root_2 = SecE4_curve_center + 2
			SecE4_curve_root_2 = self.blueprint.segment_dict['E2'].bp_data[-1][0]

			#
			# Add lines for actual cst string
			CST_line_SecE4 = "Angle CA %d CA %d CA %d HARMONIC %0.2f %0.2f" \
			%(SecE4_curve_root_1,SecE4_curve_center,SecE4_curve_root_2,mean_angle_value,std_dev)
			#
			CST_RETURN_LIST.append(CST_line_SecE4)

			# Now for E5:
			relative_SecE5_curve_center = -1*(relative_MainE4_curve_center + 1) - 5 # -4
			SecE5_curve_center = self.blueprint.segment_dict['E3'].bp_data[relative_SecE5_curve_center][0]
			#SecE5_curve_root_1 = SecE5_curve_center - 2
			SecE5_curve_root_1 = self.blueprint.segment_dict['E3'].bp_data[0][0]
			SecE5_curve_root_2 = SecE5_curve_center + 2


			CST_line_SecE5 = "Angle CA %d CA %d CA %d HARMONIC %0.2f %0.2f" \
			%(SecE5_curve_root_1,SecE5_curve_center,SecE5_curve_root_2,mean_angle_value,std_dev)
			#
			CST_RETURN_LIST.append(CST_line_SecE5)

			'''
			# DIHEDRAL:
			second_to_last_E6 = self.blueprint.bp_data[-2][0]
			pre_Mainbulge = self.sheet_data['main_bulge_E3_pos'] - 1
			second_res = 2
			first_res_E5 = self.blueprint.segment_dict['E3'].bp_data[0][0]

			mean_dihedral_value =  math.radians(80)

			CST_line_DIH = "Dihedral CA %d CA %d CA %d CA %d HARMONIC %0.2f %0.2f" \
			%(second_to_last_E6,pre_Mainbulge,second_res,first_res_E5,mean_dihedral_value,std_dev)

			if (self.long_arm_l == 4):
				relative_SecE5_curve_center = E4_len - relative_SecE4_curve_center -2 # keep at the same level?? Yes, works better.
				SecE5_curve_center = self.blueprint.segment_dict['E3'].bp_data[relative_SecE5_curve_center][0]
				SecE5_curve_root_1 = SecE5_curve_center - 2
				SecE5_curve_root_2 = SecE5_curve_center + 2
				CST_line_SecE5 = "Angle CA %d CA %d CA %d HARMONIC %0.2f %0.2f" \
				%(SecE5_curve_root_1,SecE5_curve_center,SecE5_curve_root_2,mean_angle_value,std_dev)

				CST_RETURN_LIST.append(CST_line_SecE5)
		'''
		self.sheet_data["CST_lines"] = CST_RETURN_LIST
		return CST_RETURN_LIST

	def produce_bp(self):
		# Compute length of E3:
		E3_l = self.base_width + 1 + 2*self.long_arm_l
		if self.Second_bulge_E3 : E3_l = E3_l + 1
		if self.ExtendedE4 : E3_l = E3_l + 2
		self.sheet_data["E3_l"] = E3_l
		# Compute length of E4:
		E4_l = self.base_width + 2*self.long_arm_l + 1
		if self.ExtendedE4 : E4_l = E4_l + 2
		self.sheet_data["E4_l"] = E4_l
		# Compute length of E5:
		E5_l = self.base_width + 2*self.long_arm_l + 1 + 2*self.short_arm_l
		self.sheet_data["E5_l"] = E5_l
		# Compute length of E6:
		E6_l = 2 + 2*self.short_arm_l +  self.base_width
		if self.ExtendedE6:
			E6_l = E6_l + 2
		self.sheet_data["E6_l"] = E6_l
		# Strand pairings:
		SSPAIR_3_4 = -1
		if self.ExtendedE4 :
			SSPAIR_4_5 = -2
		else:
			SSPAIR_4_5 = -2*self.short_arm_l

		SSPAIR_5_6 = E5_l - E6_l + 1
		self.pairings = ["1-2.A.%d"%SSPAIR_3_4,"2-3.A.%d"%SSPAIR_4_5,"3-4.A.%d"%SSPAIR_5_6]
		self.sheet_data["pairings"] = { "SSPAIR_3_4":SSPAIR_3_4, "SSPAIR_4_5":SSPAIR_4_5, "SSPAIR_5_6":SSPAIR_5_6 }
		self.sheet_data["pairing_lines"] = self.pairings
		# Segment creation
		L1 = Segment('L1','L',[[1,'A','LX','.']])
		segments = [L1]

		E1 = Segment('E1','E',[[i+1,'A','EB','.'] for i in range(E3_l)])
		# E3 bulge placemanets
		if self.ExtendedE4:
			main_bulge_pos_from_C = -1*(self.base_width)-3
		else:
			main_bulge_pos_from_C = -1*(self.base_width)-1

		E1.bp_data[main_bulge_pos_from_C][2] = 'EA'

		if self.Second_bulge_E3 :
			second_bulge_position_from_C = main_bulge_pos_from_C - ( 2*self.Second_b_place + 3 )
			E1.bp_data[second_bulge_position_from_C][2] = 'EA'

		segments += [E1]
		L2 = Segment('L2','L',[[1,'A','LX','.'],[2,'A','LX','.']])
		segments += [L2]
		# E4
		E2 = Segment('E2','E',[[i+1,'A','EB','.'] for i in range(E4_l)])
		segments += [E2]
		L3 = Segment('L3','L',[[1,'A','LX','.'],[2,'A','LX','.']])
		segments += [L3]
		# E5
		E3 = Segment('E3','E',[[i+1,'A','EB','.'] for i in range(E5_l)])
		segments += [E3]
		if self.short_arm_l == 2:
			L4 = Segment('L4','L',[[1,'A','LG','.'],[2,'A','LG','.']])
		else:
			L4 = Segment('L4','L',[[1,'A','LX','.'],[2,'A','LX','.']])
		segments += [L4]
		# E6
		E4 = Segment('E4','E',[[i+1,'A','EB','.'] for i in range(E6_l)])
		# E6 bulge placemanets
		E6_bulge_INDEX_from_N = 2*(self.short_arm_l)
		E4.bp_data[E6_bulge_INDEX_from_N][2] = 'EA'
		segments += [E4]
		L5 = Segment('L5','L',[[1,'A','LX','.']])
		segments += [L5]

		bp = Blueprint(segments=segments)
		bp.reindex_blueprint()
		self.sheet_data["main_bulge_E6_pos"] = bp.segment_dict['E4'].bp_data[E6_bulge_INDEX_from_N][0]
		self.sheet_data["main_bulge_E3_pos"] = bp.segment_dict['E1'].bp_data[main_bulge_pos_from_C][0]

		def aux_get_MainBulge_hinges(bp):
			if self.ExtendedE4:
				relative_MainE4_curve_center = self.base_width + 2
			else:
				relative_MainE4_curve_center = self.base_width
			E4_len = len(bp.segment_dict['E2'].bp_data)
			relative_MainE5_curve_center = E4_len - relative_MainE4_curve_center - 1
			return (relative_MainE4_curve_center,relative_MainE5_curve_center)

		def aux_get_SecBulge_hinges(bp):
			E4_len = len(bp.segment_dict['E2'].bp_data)
			positions = aux_get_MainBulge_hinges(bp)
			relative_SecE4_curve_center = positions[0] + 2*self.Second_b_place + 2
			relative_SecE5_curve_center = E4_len - relative_SecE4_curve_center - 1
			return (relative_SecE4_curve_center,relative_SecE5_curve_center)

		if self.Second_bulge_E3 :
			self.sheet_data["sec_bulge_E3_pos"] =  bp.segment_dict['E1'].bp_data[second_bulge_position_from_C][0]
			if self.E3_SecBulgeCurve == 3:
				hinges_pos_sec = aux_get_SecBulge_hinges(bp)
				gly_pos = hinges_pos_sec[0]
				#bp.segment_dict['E2'].bp_data[gly_pos][1] = 'G'
		if self.E3_MainBulgeCurve == 3:
			hinges_pos_main = aux_get_MainBulge_hinges(bp)
			gly_pos = hinges_pos_main[0]
			#bp.segment_dict['E2'].bp_data[gly_pos][1] = 'G'
		return bp

	def get_MainBulge_hinges(self):
		bp = self.blueprint
		if self.ExtendedE4:
			relative_MainE4_curve_center = self.base_width + 2
		else:
			relative_MainE4_curve_center = self.base_width
		E4_len = len(bp.segment_dict['E2'].bp_data)
		relative_MainE5_curve_center = E4_len - relative_MainE4_curve_center - 1
		return (relative_MainE4_curve_center,relative_MainE5_curve_center)

	def get_SecBulge_hinges(self):
		bp = self.blueprint
		E4_len = len(bp.segment_dict['E2'].bp_data)
		positions = self.get_MainBulge_hinges()
		relative_SecE4_curve_center = positions[0] + 2*self.Second_b_place + 2
		relative_SecE5_curve_center = E4_len - relative_SecE4_curve_center - 1
		return (relative_SecE4_curve_center,relative_SecE5_curve_center)

	def get_CurvedLongArm_hinges(self):
		bp = self.blueprint
		if self.ExtendedE4:
			relative_MainE4_curve_center = self.base_width + 2
		else:
			relative_MainE4_curve_center = self.base_width
		#main_hinge_rel_pos = self.get_MainBulge_hinges()
		relative_SecE4_curve_center = relative_MainE4_curve_center + 4
		#SecE4_curve_center = self.blueprint.segment_dict['E2'].bp_data[relative_SecE4_curve_center-1][0]
		return relative_SecE4_curve_center

	def write_blueprint(self,fname=None):
		if not fname:
			fname = './sheet.bp'
		sspairs_str = "SSPAIR "+";".join(self.pairings)
		self.blueprint.dump_blueprint(fname,[sspairs_str])

	def write_cstfile(self,fname=None):
		if not fname:
			fname = './angle.csts'
		outfile = open(fname,'w')
		for line in self.csts:
			outfile.write(line+'\n')
		outfile.close()

	def get_Tomponents_flags_dict(self):
		tomponents_dict = {'E3_l':self.sheet_data["E3_l"]}
		tomponents_dict['E4_l'] = self.sheet_data["E4_l"]
		tomponents_dict['E5_l'] = self.sheet_data["E5_l"]
		tomponents_dict['E6_l'] = self.sheet_data["E6_l"]
		E6_bulge_relative = self.sheet_data["main_bulge_E6_pos"] - self.blueprint.segment_dict['E4'].bp_data[0][0] + 1
		tomponents_dict['BulgeE6_str'] = '%d'%E6_bulge_relative
		E3_n = self.blueprint.segment_dict['E1'].bp_data[0][0]
		main_E3_bulge_relative = self.sheet_data["main_bulge_E3_pos"] - E3_n + 1
		if self.Second_bulge_E3:
			sec_E3_bulge_relative = self.sheet_data["sec_bulge_E3_pos"] - E3_n + 1
			tomponents_dict['BulgeE3_str'] = "%d,%d"%(main_E3_bulge_relative,sec_E3_bulge_relative)
		else:
			tomponents_dict['BulgeE3_str'] = "%d"%main_E3_bulge_relative
		tomponents_dict["RS_1"] = "%d"%(-1*self.sheet_data["pairings"]["SSPAIR_3_4"])
		#tomponents_dict["RS_1"] = "0"
		tomponents_dict["RS_2"] = "%d"%(-1*self.sheet_data["pairings"]["SSPAIR_4_5"])
		tomponents_dict["RS_3"] = "%d"%self.sheet_data["pairings"]["SSPAIR_5_6"]

		return tomponents_dict

	def write_Tomponents_flags(self,fname=None):
		if not fname:
			fname = './tomponents.vals'
		outfile = open(fname,'w')
		outfile.write("-parser:script_vars E3_l=%d \n"%self.sheet_data["E3_l"])
		outfile.write("-parser:script_vars E4_l=%d \n"%self.sheet_data["E4_l"])
		outfile.write("-parser:script_vars E5_l=%d \n"%self.sheet_data["E5_l"])
		outfile.write("-parser:script_vars E6_l=%d \n"%self.sheet_data["E6_l"])
		E6_bulge_relative = self.sheet_data["main_bulge_E6_pos"] - self.blueprint.segment_dict['E4'].bp_data[0][0] + 1
		outfile.write("-parser:script_vars BulgeE6_str=%d \n"%E6_bulge_relative)
		E3_n = self.blueprint.segment_dict['E1'].bp_data[0][0]
		main_E3_bulge_relative = self.sheet_data["main_bulge_E3_pos"] - E3_n + 1
		if self.Second_bulge_E3:
			sec_E3_bulge_relative = self.sheet_data["sec_bulge_E3_pos"] - E3_n + 1
			outfile.write("-parser:script_vars BulgeE3_str=%d,%d \n"%(main_E3_bulge_relative,sec_E3_bulge_relative) )
		else:

			outfile.write("-parser:script_vars BulgeE3_str=%d \n"%main_E3_bulge_relative )
		outfile.write("-parser:script_vars RS_1=%d \n"%self.sheet_data["pairings"]["SSPAIR_3_4"])
		outfile.write("-parser:script_vars RS_2=%d \n"%(-1*self.sheet_data["pairings"]["SSPAIR_4_5"]))
		outfile.write("-parser:script_vars RS_3=%d \n"%self.sheet_data["pairings"]["SSPAIR_5_6"])
		outfile.close()

	def dump_sheet_data_human_readable(self,fname=None):
		if not fname:
			fname = './sheet_data.txt'
		outfile = open(fname,'w')
		for i in self.sheet_data.keys():
			outfile.write("%s %s\n"%(i,self.sheet_data[i]))
		outfile.close()

def Get_SheetE6bulge_E3N_dist(sheet_pose,sheet):
	sheet_bp = Blueprint(data=[ [ x for x in i ] for i in sheet.blueprint.bp_data ] )
	sheet_bp.reindex_blueprint()
	res1_pos = sheet_bp.segment_dict['E1'].bp_data[0][0]
	res2_pos = [i[0]+1 for i in sheet_bp.segment_dict["E4"].bp_data if i[2]=='EA'][0]
	CA_1 = sheet_pose.residue(res1_pos).xyz('CA')
	CA_2 = sheet_pose.residue(res2_pos).xyz('CA')
	CA_CAv = CA_2 - CA_1
	return CA_CAv.length()

def Get_SheetBase_longArm_angle(sheet_pose,sheet):
	sheet_bp = Blueprint(data=[ [ x for x in i ] for i in sheet.blueprint.bp_data ] )
	sheet_bp.reindex_blueprint()
	# Get atoms and their coordinates for plane formed by tip of long arm
	lArm_E3N = sheet_bp.segment_dict['E1'].bp_data[0][0]
	lArm_E3N_2 = lArm_E3N+1
	lArm_E5N = sheet_bp.segment_dict['E3'].bp_data[0][0]
	lArm_E3N_CA = sheet_pose.residue(lArm_E3N).xyz('CA')
	lArm_E3N_2_CA = sheet_pose.residue(lArm_E3N_2).xyz('CA')
	lArm_E5N_CA = sheet_pose.residue(lArm_E5N).xyz('CA')
	cys_v_LA = lArm_E3N_2_CA - lArm_E3N_CA
	trans_v_LA = lArm_E5N_CA - lArm_E3N_CA
	# Get an orthogonal vector to the long arm plane
	orth_LA_v = cys_v_LA.cross_product(trans_v_LA)
	# Get atoms and their coordinates for plane formed by the sheet base
	base_E6_C_1 = sheet_bp.segment_dict['E4'].bp_data[-2][0]
	base_E6_bulge = sheet.sheet_data["main_bulge_E6_pos"]
	base_E3_bulge_1 = sheet.sheet_data["main_bulge_E3_pos"]+1
	base_E6_C_1_CA = sheet_pose.residue(base_E6_C_1).xyz('CA')
	base_E6_bulge_CA = sheet_pose.residue(base_E6_bulge).xyz('CA')
	base_E3_bulge_1_CA = sheet_pose.residue(base_E3_bulge_1).xyz('CA')
	cys_v_base = base_E6_bulge_CA - base_E6_C_1_CA
	trans_v_base = base_E3_bulge_1_CA - base_E6_C_1_CA
	# Get an orthogonal vector to the long arm plane
	orth_base_v = cys_v_base.cross_product(trans_v_base)
	# Finally, calculate the angles between the base plane and the long arm tip plane:
	N_orth_LA_v = orth_LA_v.normalize_or_zero()
	N_orth_base_v = orth_base_v.normalize_or_zero()
	angle = math.degrees( math.acos( rosetta.numeric.sin_cos_range(N_orth_base_v.dot_product(N_orth_LA_v)) ) )
	return angle

def Get_longArm_shootingAngle(sheet_pose,sheet):
	sheet_bp = Blueprint(data=[ [ x for x in i ] for i in sheet.blueprint.bp_data ] )
	sheet_bp.reindex_blueprint()
	# Get atoms and their coordinates for plane formed by the sheet base
	base_E6_C_1 = sheet_bp.segment_dict['E4'].bp_data[-2][0]
	base_E6_bulge = sheet.sheet_data["main_bulge_E6_pos"]
	base_E3_bulge_1 = sheet.sheet_data["main_bulge_E3_pos"]+1
	base_E6_C_1_CA = sheet_pose.residue(base_E6_C_1).xyz('CA')
	base_E6_bulge_CA = sheet_pose.residue(base_E6_bulge).xyz('CA')
	base_E3_bulge_1_CA = sheet_pose.residue(base_E3_bulge_1).xyz('CA')
	cys_v_base = base_E6_bulge_CA - base_E6_C_1_CA
	trans_v_base = base_E3_bulge_1_CA - base_E6_C_1_CA
	# Get an orthogonal vector to the base plane
	orth_base_v = cys_v_base.cross_product(trans_v_base)
	# Get vector of long arm shooting out:
	lArm_E3N = sheet_bp.segment_dict['E1'].bp_data[0][0]
	lArm_E3N_2 = lArm_E3N+1
	lArm_E3N_CA = sheet_pose.residue(lArm_E3N).xyz('CA')
	lArm_E3N_2_CA = sheet_pose.residue(lArm_E3N_2).xyz('CA')
	cys_v_LA = lArm_E3N_CA - lArm_E3N_2_CA
	#Angle:
	N_cys_v_LA = cys_v_LA.normalize_or_zero()
	N_orth_base_v = orth_base_v.normalize_or_zero()
	angle = math.degrees( math.acos( rosetta.numeric.sin_cos_range(N_cys_v_LA.dot_product(N_orth_base_v)) ) )
	return 180-angle

def Get_longArm_protrusionDist(sheet_pose,sheet):
	sheet_bp = Blueprint(data=[ [ x for x in i ] for i in sheet.blueprint.bp_data ] )
	sheet_bp.reindex_blueprint()
	base_E3_bulge_1 = sheet.sheet_data["main_bulge_E3_pos"]+1
	base_E6_ref = sheet.sheet_data["main_bulge_E6_pos"]+sheet.sheet_data["base_width"]
	base_E3_lArm_N = sheet_bp.segment_dict['E1'].bp_data[0][0]
	base_E6_CA = sheet_pose.residue(base_E6_ref).xyz('CA')
	base_E3bulge_CA = sheet_pose.residue(base_E3_bulge_1).xyz('CA')
	lArmN_CA = sheet_pose.residue(base_E3_lArm_N).xyz('CA')
	a = lArmN_CA - base_E3bulge_CA
	b = base_E6_CA - base_E3bulge_CA
	projection = a.dot(b)/b.length()
	protrusion = projection - b.length()
	return protrusion

def GetSheetType(sheet_obj,db):
	main_dict = sheet_obj.saved_input
	# Generate the sheet type DB to confront with:
	main_db = {}
	for i in range(1,69):
		#fname = '/home/basantab/NTF2_project/20170706_CreateRingFromSheet/AutomaticRingSelectionTest/SheetTypeDB/%04d_sheet.dict'%i
		fname = '%s/all_sheet_dicts/%04d_sheet.dict'%(db,i)
		main_db[i] = json.load(open(fname,'r'))
	for key in main_db.keys():
		if all( [ (main_dict[i] ==  main_db[key][i]) for i in main_db[key].keys()] ):
			return key
		else:
			continue
	raise ValueError('NO SHEET TYPE FOUND')

class RingConnection():

	Connection_types = ['BA','GBA','GB','ClassicDirect','BulgeAndB','BBGB']

	def generate_bp(self):

		if self.connection_type == 'BulgeAndB':
			# Make 2 first residues remodel
			self.ring_bp.bp_data[0][-1] = 'R'
			self.ring_bp.bp_data[1][-1] = 'R'
			self.ring_bp.bp_data[2][-1] = 'R'
			# Add residues on N term:
			self.ring_bp.bp_data[0][-2] = 'EB'
			self.pairings.append('3-4.A.0')
			self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LA','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LG','R'])
			for i in range(self.h_len):
				self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])

		if self.connection_type == 'BA': # This is a "Classic Bulge"
			# Make 2 first residues remodel
			self.ring_bp.bp_data[0][-1] = 'R'
			self.ring_bp.bp_data[1][-1] = 'R'
			self.ring_bp.bp_data[2][-1] = 'R'
			# Add residues on N term:
			self.ring_bp.bp_data[0][-2] = 'LA'
			self.pairings.append('3-4.A.-1')
			self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
			for i in range(self.h_len):
				self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])
		if self.connection_type == 'ClassicDirect':
			# Make 3 first residues remodel
			self.ring_bp.bp_data[0][-1] = 'R'
			self.ring_bp.bp_data[1][-1] = 'R'
			self.ring_bp.bp_data[2][-1] = 'R'
			self.ring_bp.bp_data[0][-2] = 'LB'
			self.pairings.append('3-4.A.-1')
			#self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
			# Add residues on N term:
			#self.ring_bp.bp_data[0][-2] = 'HA'
			for i in range(self.h_len):
				self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])
		if self.connection_type == 'ClassicBulge':
			# Make 3 first residues remodel
			self.ring_bp.bp_data[0][-1] = 'R'
			self.ring_bp.bp_data[1][-1] = 'R'
			self.ring_bp.bp_data[2][-1] = 'R'
			#self.ring_bp.bp_data[3][-1] = 'R'
			#self.ring_bp.bp_data[4][-1] = 'R'
			# Add residues on N term:
			self.ring_bp.bp_data[1][-2] = 'EA'
			self.ring_bp.bp_data[0][-2] = 'LB'
			self.pairings.append('3-4.A.-1')
			for i in range(self.h_len):
				self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])
		if self.connection_type == 'GBA':
			self.ring_bp.bp_data[0][-1] = 'R'
			self.ring_bp.bp_data[1][-1] = 'R'
			self.ring_bp.bp_data[2][-1] = 'R'
			self.ring_bp.bp_data[0][-2] = 'LA'
			self.pairings.append('3-4.A.-1')
			self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LG','R'])
			for i in range(self.h_len):
				self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])
		if self.connection_type == 'GB':
			self.ring_bp.bp_data[0][-1] = 'R'
			self.ring_bp.bp_data[1][-1] = 'R'
			self.ring_bp.bp_data[2][-1] = 'R'
			self.ring_bp.bp_data[0][-2] = 'LB'
			self.pairings.append('3-4.A.-1')
			self.ring_bp.bp_data.insert(0,[0,'A','LG','R'])
			for i in range(self.h_len):
				self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])
		if self.connection_type == 'BBGB':
			self.ring_bp.bp_data[0][-1] = 'R'
			self.ring_bp.bp_data[1][-1] = 'R'
			#self.ring_bp.bp_data[1][-2] = 'B'
			self.ring_bp.bp_data[0][-2] = 'LB'
			self.pairings.append('3-4.A.-1')
			self.ring_bp.bp_data.insert(0,[0,'A','LG','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
			for i in range(self.h_len):
				self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])
		if self.connection_type == 'GAA':
			self.ring_bp.bp_data[0][-1] = 'R'
			self.ring_bp.bp_data[1][-1] = 'R'
			self.ring_bp.bp_data[1][-2] = 'EB'
			self.ring_bp.bp_data[0][-2] = 'LA'
			self.pairings.append('3-4.A.-1')
			self.ring_bp.bp_data.insert(0,[0,'A','LA','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LG','R'])
			for i in range(self.h_len):
				self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])
		if self.connection_type == 'BABAB':
			self.ring_bp.bp_data[0][-1] = 'R'
			self.ring_bp.bp_data[1][-1] = 'R'
			self.ring_bp.bp_data[1][-2] = 'EB'
			self.ring_bp.bp_data[0][-2] = 'LB'
			self.pairings.append('3-4.A.-1')
			self.ring_bp.bp_data.insert(0,[0,'A','LA','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LA','R'])
			self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
			for i in range(self.h_len):
				self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])

		############## COMMMON PART ####################
		self.ring_bp.bp_data.insert(0,[0,'A','L%s'%self.loopOneABEGO,'R'])
		for i in range(self.hairpin_len):
			self.ring_bp.bp_data.insert(0,[0,'A','EB','R'])
		self.ring_bp.bp_data.insert(0,[0,'A','LX','R'])
		self.ring_bp.bp_data.insert(0,[0,'A','LX','R'])
		for i in range(self.hairpin_len):
			self.ring_bp.bp_data.insert(0,[0,'A','EB','R'])
		self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
		self.ring_bp.bp_data.insert(0,[0,'A','LA','R'])
		self.ring_bp.bp_data.insert(0,[0,'A','LA','R'])
		self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
		self.ring_bp.bp_data.insert(0,[0,'A','LB','R'])
		for i in range(7):
			self.ring_bp.bp_data.insert(0,[0,'A','HA','R'])
		self.ring_bp.bp_data.insert(0,[0,'A','LX','R'])
		self.pairings.append('1-2.A.0')
		#self.pairings.append('3-4.A.0')
		E6_E1_RegShift = -1*(self.sheet_obj.sheet_data["short_arm_l"]*2 + 1)
		hairpin_pairing_string = '1-6.P.%d'%E6_E1_RegShift
		self.pairings.append(hairpin_pairing_string)
		#self.ring_bp.freeze_all()
		#self.ring_bp.reindex_blueprint()

	def populate_ring_dict(self):
		self.ring_dict ={
		'h_len':self.h_len, \
		'hairpin_len':self.hairpin_len, \
		'connection_type':self.connection_type,\
		'loopOneABEGO':self.loopOneABEGO, \
		'sheet_dict':self.sheet_obj.sheet_data, \
		'pairings':self.pairings, \
		'sheet_type':GetSheetType(self.sheet_obj,self.db) \
		}

	ring_attributes = ['h_len','hairpin_len','connection_type','loopOneABEGO']

	def __init__(self, sheet_obj, sheet_pose = None , h_len = None, loopOneABEGO = 'E', hairpin_len = 4, connection_type='GB',\
			db = None ,bias_TP=False, bias_CH_not_TP=False):
		self.allowed_connections = self.Connection_types
		self.db = db
		self.h_len = h_len
		self.hairpin_len = hairpin_len
		if connection_type not in self.allowed_connections:
			raise ValueError('Connection type for ring not recognized, allowed values are: ABG, BAA, BAG, BBB, GAA, ClassicBulge, ClassicDirect, BABAB and BulgeAndB ')
		self.connection_type = connection_type
		self.loopOneABEGO = loopOneABEGO
		self.sheet_pose = sheet_pose
		self.sheet_obj = sheet_obj
		self.bias_TP = bias_TP
		self.bias_CH_not_TP = bias_CH_not_TP
		if self.bias_TP and self.bias_CH_not_TP: raise ValueError('You can\'t bias for C-terminal helix and TP at the same time')
		self.ring_bp = Blueprint(data=[ [ x for x in i ] for i in sheet_obj.blueprint.bp_data ] )
		self.ring_bp.freeze_all()
		self.ring_bp.reindex_blueprint()
		self.pairings = []
		inherited_pairings = [ "4-5.A.%s"%(i.split('.')[-1]) for i in self.sheet_obj.sheet_data['pairing_lines'] if "2-3.A" in i ]
		inherited_pairings += [ "5-6.A.%s"%(i.split('.')[-1]) for i in self.sheet_obj.sheet_data['pairing_lines'] if "3-4.A" in i ]
		self.pairings += inherited_pairings
		self.generate_bp()
		self.ring_dict = {}
		self.populate_ring_dict()

	def write_blueprint(self,fname=None):
		if not fname:
			fname = './ring.bp'
		sspairs_str = "SSPAIR "+";".join(self.pairings)
		self.ring_bp.dump_blueprint(fname,[sspairs_str])

	def get_min_fix_res(self):
		'''
		Get the residues to fix during minimization step.
		'''
		last = len(self.ring_bp.bp_data)
		first = 1
		for n,i in enumerate(self.ring_bp.bp_data):
			if i[-1] !='R':
				first = n+1
				break
		return (first,last)
	def dump_all_csts(self,fname=None):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.ring_bp.bp_data ] )
		bp.reindex_blueprint()

		cst_string = ''
		if self.bias_TP:
			E4c = bp.segment_dict['E4'].bp_data[-1][0]
			E1c = bp.segment_dict['E1'].bp_data[-1][0]
			cst_string += "AtomPair CA %i CA %i FLAT_HARMONIC 14.5 3.0 1.4\n"%(E1c,E4c)

		if self.bias_CH_not_TP:
			E6c = bp.segment_dict['E6'].bp_data[-1][0]
			E5n = bp.segment_dict['E5'].bp_data[0][0]
			cst_string += "AtomPair CA %i CA %i FLAT_HARMONIC 16.5 3.0 1.9\n"%(E6c,E5n)

		if self.connection_type=='ClassicBulge':
			hbond1_don = bp.segment_dict['E3'].bp_data[0][0]
			hbond1_acc = bp.segment_dict['E4'].bp_data[-2][0]
			cst_string += CircularHBondConstraints(hbond1_don,hbond1_acc)

		cst_string += PerfectHelixCst(self.ring_bp,1)
		cst_string += PerfectHelixCst(self.ring_bp,2)

		hbond1_don = bp.segment_dict['L4'].bp_data[0][0]
		hbond1_acc = bp.segment_dict['L2'].bp_data[-1][0]
		cst_string += CircularHBondConstraints(hbond1_don,hbond1_acc)

		hbond2_don = bp.segment_dict['H2'].bp_data[0][0]
		hbond2_acc = bp.segment_dict['L2'].bp_data[-3][0]
		cst_string += CircularHBondConstraints(hbond2_don,hbond2_acc)

		hbond3_don = bp.segment_dict['L2'].bp_data[-1][0]
		hbond3_acc = bp.segment_dict['L2'].bp_data[-4][0]
		cst_string += CircularHBondConstraints(hbond3_don,hbond3_acc)

		hbond4_don = bp.segment_dict['L2'].bp_data[-4][0]
		hbond4_acc_pos = 2*(self.sheet_obj.sheet_data["short_arm_l"]) - 1
		hbond4_acc = bp.segment_dict['E6'].bp_data[hbond4_acc_pos][0]
		cst_string += CircularHBondConstraints(hbond4_don,hbond4_acc)

		hbond5_acc = bp.segment_dict['H1'].bp_data[-1][0]
		hbond5_don_pos = 2*(self.sheet_obj.sheet_data["short_arm_l"]) - 1
		hbond5_don = bp.segment_dict['E6'].bp_data[hbond5_don_pos][0]
		cst_string += CircularHBondConstraints(hbond5_don,hbond5_acc)
		if not fname:
			fname = './ring.csts'
		handle = open(fname,'w')
		handle.write(cst_string)
		handle.close()
	def dump_Hbonds_csts(self,fname=None):
		bp = Blueprint( data=[ [ x for x in i ] for i in self.ring_bp.bp_data ] )
		bp.reindex_blueprint()
		cst_string = ''

		hbond1_don = bp.segment_dict['L4'].bp_data[0][0]
		hbond1_acc = bp.segment_dict['L2'].bp_data[-1][0]
		cst_string += CircularHBondConstraints(hbond1_don,hbond1_acc)

		hbond2_don = bp.segment_dict['H2'].bp_data[0][0]
		hbond2_acc = bp.segment_dict['L2'].bp_data[-3][0]
		cst_string += CircularHBondConstraints(hbond2_don,hbond2_acc)

		hbond3_don = bp.segment_dict['L2'].bp_data[-1][0]
		hbond3_acc = bp.segment_dict['L2'].bp_data[-4][0]
		cst_string += CircularHBondConstraints(hbond3_don,hbond3_acc)

		hbond4_don = bp.segment_dict['L2'].bp_data[-4][0]
		hbond4_acc_pos = 2*(self.sheet_obj.sheet_data["short_arm_l"]) - 1
		hbond4_acc = bp.segment_dict['E6'].bp_data[hbond4_acc_pos][0]
		cst_string += CircularHBondConstraints(hbond4_don,hbond4_acc)

		hbond5_acc = bp.segment_dict['H1'].bp_data[-1][0]
		hbond5_don_pos = 2*(self.sheet_obj.sheet_data["short_arm_l"]) - 1
		hbond5_don = bp.segment_dict['E6'].bp_data[hbond4_acc_pos][0]
		cst_string += CircularHBondConstraints(hbond5_don,hbond5_acc)
		if not fname:
			fname = './ring_hbond.csts'
		handle = open(fname,'w')
		handle.write(cst_string)
		handle.close()

	def get_noPro_posL2(self):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.ring_bp.bp_data ] )
		bp.reindex_blueprint()
		l2 = bp.segment_dict['L2'].bp_data[1][0]
		l2c = bp.segment_dict['L2'].bp_data[-1][0]
		h3n = bp.segment_dict['H2'].bp_data[0][0]
		E3n = bp.segment_dict['E3'].bp_data[0][0]
		return '%d,%d,%d,%d'%(l2,h3n,l2c,E3n)

	def FrontHPHbonds(self):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.ring_bp.bp_data ] )
		bp.reindex_blueprint()
		return_list = []
		#BulgePos_rel = 2*(self.sheet_obj.sheet_data["short_arm_l"])
		#BulgePos = bp.segment_dict['E6'].bp_data[BulgePos_rel][0]
		#return_list.append(BulgePos+2)
		#return_list.append(BulgePos+4)
		HP_n = bp.segment_dict['E1'].bp_data[0][0]
		return_list.append(HP_n)
		return_list.append(HP_n+2)

		#E6 = [ i[0] for i in bp.segment_dict['E6'].bp_data ]
		#E1 = [ i[0] for i in bp.segment_dict['E1'].bp_data ]

		#return_list = E6 + E1
		#if self.hairpin_len >= 6:
		#	return_list.append(BulgePos+6)
		#	return_list.append(HP_n+4)
		return return_list

def Get_RinglongArm_protrusionDist(ring_pose,ring_obj):
	ring_bp = Blueprint(data=[ [ x for x in i ] for i in ring_obj.ring_bp.bp_data ] )
	added_aa_n = len( [i for i in ring_bp.bp_data if i[0] == 0 ] )
	base_E3_bulge_1 = ring_obj.ring_dict['sheet_dict']["main_bulge_E3_pos"]+1 + added_aa_n
	base_E6_ref = ring_obj.ring_dict['sheet_dict']["main_bulge_E6_pos"]+ring_obj.ring_dict['sheet_dict']["base_width"] + added_aa_n
	base_E3_lArm_N = ring_bp.segment_dict['E3'].bp_data[0][0] + added_aa_n
	base_E6_CA = ring_pose.residue(base_E6_ref).xyz('CA')
	base_E3bulge_CA = ring_pose.residue(base_E3_bulge_1).xyz('CA')
	lArmN_CA = ring_pose.residue(base_E3_lArm_N).xyz('CA')
	a = lArmN_CA - base_E3bulge_CA
	b = base_E6_CA - base_E3bulge_CA
	projection = a.dot(b)/b.length()
	protrusion = projection - b.length()
	return protrusion

def CreateSheetObjFromDict(sheet_data_dict):
	'''
	Creates sheet from input dictionary, that dictionaty can be stored sheet_data.
	'''
	dummy_dict = {}
	for i in NTF2_sheet.sheet_atributes_Dict.keys():
		dummy_dict[i] = sheet_data_dict[i]

	return NTF2_sheet(**dummy_dict)

def CreateRingObjFromDict(ring_data_dict,db=None,bias_TP=False,bias_CH_not_TP=False):
	'''
	Creates a ring from input dictionary, that dictionaty can be stored ring_dict.
	'''
	dummy_dict = {}
	for i in RingConnection.ring_attributes:
		dummy_dict[i] = ring_data_dict[i]
	sheet_obj = CreateSheetObjFromDict(ring_data_dict['sheet_dict'])

	return RingConnection(sheet_obj,**dummy_dict,db=db,bias_TP=bias_TP,bias_CH_not_TP=bias_CH_not_TP)

def CreateAllPossibleSheetDicts():
	all_sheet_types = []
	for var in NTF2_sheet.sheet_atributes_Dict['base_width']:
		for var1 in NTF2_sheet.sheet_atributes_Dict['long_arm_l']:
			for var2 in NTF2_sheet.sheet_atributes_Dict['short_arm_l']:
				if var2==2:
					extendedE4=True
				else:
					extendedE4=False
				for var3 in NTF2_sheet.sheet_atributes_Dict['E3_MainBulgeCurve']:
					for var4 in NTF2_sheet.sheet_atributes_Dict['ExtendedE6']:
						for var5 in NTF2_sheet.sheet_atributes_Dict['Second_bulge_E3']:
							if var5:
								for var6 in NTF2_sheet.sheet_atributes_Dict['Second_b_place']:
									for var7 in NTF2_sheet.sheet_atributes_Dict['E3_SecBulgeCurve']:
										current = {
											"base_width" : var, \
											"long_arm_l" : var1, \
											"short_arm_l" : var2, \
											"Second_bulge_E3" : var5, \
											"Second_b_place" : var6, \
											"E3_MainBulgeCurve" : var3, \
											"E3_SecBulgeCurve" : var7, \
											"ExtendedE6" : var4, \
											"ExtendedE4" : extendedE4
										}
										sheet =  NTF2_sheet(**current, validate_sheet=False)
										try:
											sheet.validate_sheet()
										except ValueError as e:
											continue
										if sheet.sheet_data["Validated"]:
											all_sheet_types.append(current)
							else:
								for var8 in NTF2_sheet.sheet_atributes_Dict["CurvedLongArm"]:
									current = {
										"base_width" : var, \
										"long_arm_l" : var1, \
										"short_arm_l" : var2, \
										"Second_bulge_E3" : var5, \
										"E3_MainBulgeCurve" : var3, \
										"ExtendedE6" : var4, \
										"ExtendedE4" : extendedE4, \
										"CurvedLongArm" : var8 \
									}
									sheet =  NTF2_sheet(**current, validate_sheet=False)
									try:
										sheet.validate_sheet()
									except ValueError as e:
										continue
									if sheet.sheet_data["Validated"]:
										all_sheet_types.append(current)
	return [ i for i in enumerate(all_sheet_types,start=1) ]

def ExtendableRingHP(sheet_obj):
	extedable_base3_sheet = ( sheet_obj.sheet_data['base_width'] == 3 ) and sheet_obj.sheet_data['ExtendedE6']
	extendable_HP = extedable_base3_sheet or ( sheet_obj.sheet_data['base_width'] == 5 )
	return extendable_HP

def ParseConnections(sheet_pose, sheet_obj):
	protrusion = Get_longArm_protrusionDist(sheet_pose,sheet_obj)
	if protrusion < -7:
		return ['BulgeAndB','BBGB','GBA']
	elif protrusion > 2 or ( sheet_obj.sheet_data['Second_bulge_E3'] and protrusion > -1) : # With a smaller E3, now protrussion is shifter towards smaller values, ~2A.
		return ['BA','ClassicDirect']
	else:
		return ['BA','GBA','GB','ClassicDirect','BulgeAndB','BBGB']

def AutomaticRingRecomendation(sheet_pose, sheet_obj,db=None):

	ResultDB_fname = '%s/FinalTable.tab'%db
	ResultDB = pnd.read_csv(ResultDB_fname,delim_whitespace=True,index_col=0)
	type_list = list(set([ i for i in ResultDB.index]))
	ConnectionStatsDict = {}
	for Sheet_type in type_list:
		ConnectionStatsDict[Sheet_type] = [(i,sum(ResultDB.get_value(Sheet_type,col=i))) for i in list(ResultDB.columns)[1:]]
	sheet_type = GetSheetType(sheet_obj,db)

	#recomended_connections = [ i[0] for i in sorted(ConnectionStatsDict[sheet_type],key=lambda x: x[1],reverse=True) if i[1]!=0 ]
	recomended_connections = ['BABAB','BulgeAndB','ABG','ClassicDirect','BBB','ClassicBulge']

	protrusion = Get_longArm_protrusionDist(sheet_pose,sheet_obj)
	bulge6_dist = Get_SheetE6bulge_E3N_dist(sheet_pose,sheet_obj)
	if 'BAG' in recomended_connections: recomended_connections.remove('BAG')
	#if 'ClassicBulge' in recomended_connections: recomended_connections.remove('ClassicBulge')
	if 'GBBAA' in recomended_connections: recomended_connections.remove('GBBAA')
	if protrusion < -5:
		# These are short connections that won't work well with highly negative protrusion.
		if 'ABG' in recomended_connections: recomended_connections.remove('ABG')
		if 'BAA' in recomended_connections: recomended_connections.remove('BAA')
		if 'BBB' in recomended_connections: recomended_connections.remove('BBB')
		if 'GAA' in recomended_connections: recomended_connections.remove('GAA')
		if 'ClassicDirect' in recomended_connections: recomended_connections.remove('ClassicDirect')
	if protrusion > 5 or ( sheet_obj.sheet_data['Second_bulge_E3'] and protrusion > 1) :
		# Use only short connections for high positive protrusion:
		if 'BABAB' in recomended_connections: recomended_connections.remove('BABAB')
		if 'BulgeAndB' in recomended_connections: recomended_connections.remove('BulgeAndB')
		#if 'GBBAA' in recomended_connections: recomended_connections.remove('GBBAA')
	if bulge6_dist > 25:
		HelixChoise = [14,15]
	else:
		HelixChoise = [10,11]

	extedable_base3_sheet = ( sheet_obj.sheet_data['base_width'] == 3 ) and sheet_obj.sheet_data['ExtendedE6']
	extendable_HP = extedable_base3_sheet or ( sheet_obj.sheet_data['base_width'] == 5 )

	if extendable_HP:
		hp_len = [4,6]
	else:
		hp_len = [4]

	return (recomended_connections,HelixChoise,hp_len)

def CreateRingsFromRecomendation(recomended_connections,HelixChoises,hp_lens,sheet_obj,db):
	# sheet_pose, sheet_obj, h_len = None, loopOneABEGO = 'E', hairpin_len =
	connection_vector =  []
	for H_len in HelixChoises:
		for Conn in recomended_connections:
			for HPlen in hp_lens:
				connection_vector.append(RingConnection( sheet_obj, h_len = H_len, hairpin_len = HPlen, connection_type=Conn,db=db))
	return connection_vector

class BasicBeNTF2():

	def generate_bp(self):
		E1n_pos = self.NTF2_bp.segment_dict['E1'].bp_data[0][0]
		# Remodel all Nterm:
		for i in range(E1n_pos):
			self.NTF2_bp.bp_data[i][-1] = 'R'
		# add residues to H2 if >7
		ad_hoc_vector = []

		for i in range(self.h2_len - 7):
			ad_hoc_vector.append([0,'A','HA','R'])
		# Create H1
		for i in range(len(self.loopTwoABEGO)):
			inverse_loop_index = (-1)*(i+1)
			ad_hoc_vector.append([0,'A','L%s'%self.loopTwoABEGO[inverse_loop_index],'R'])

		for n,i in enumerate(ad_hoc_vector):
			if n == 0:
				self.NTF2_bp.bp_data[0][-1] = 'R'
				self.NTF2_bp.bp_data[0][-2] = '%s'%(i[2])
			else:
				self.NTF2_bp.bp_data.insert(0,i)

		for i in range(self.h1_len):
			self.NTF2_bp.bp_data.insert(0,[0,'A','HA','R'])
		self.NTF2_bp.bp_data.insert(0,[0,'A','LX','R'])

		# Remodel Short arm:
		# Get E6 bulge pos:
		#E6B_rel = 0
		#E6B_abs = 0
		#
		#for i,res in enumerate(self.NTF2_bp.segment_dict['E6'].bp_data):
		#	if res[-2] == 'EA':
		#		E6B_rel = i+1
		#		E6B_abs = res[0]
		#		break

		#RemodelStretchLength = 2*E6B_rel + 2
		#Startpos = E6B_abs - 2*E6B_rel - 1
		#for n,i in enumerate(self.NTF2_bp.bp_data):
		#	 if i[0] in range(Startpos,E6B_abs+1):
		#		  Uncoment to activate remodeling of short arm:
		#		 self.NTF2_bp.bp_data[n][-1] = 'R'
		#		 print('Spaceholder')


	def populate_NTF2_dict(self):
		self.NTF2_dict = {
				'ring_dict':self.ring_obj.ring_dict, \
				'h1_len':self.h1_len, \
				'h2_len':self.h2_len,\
				'Opening':self.Opening,\
				'loopTwoABEGO':self.loopTwoABEGO, \
				'protrusion':self.protrusion, \
				'has_cHelix':self.has_cHelix, \
				'c_helix_dict':self.c_helix_dict \
				}

	def __init__(self, ring_obj, ring_pose = None , h1_len = 16, loopTwoABEGO = 'GB', h2_len = 7, Opening='Classic', protrusion=0 , chelix = False):
		self.opening_types = ['Classic', 'Tropical' ]
		if h2_len < 7:
			raise ValueError('H2 less than 7 residues long is not allowed for code simplicity reasons.')
		self.h1_len = h1_len
		self.h2_len = h2_len
		if Opening not in self.opening_types:
			raise ValueError('The opening type is not recognized, allowed openings are Classic and Tropical ')
		self.Opening = Opening
		self.loopTwoABEGO = loopTwoABEGO
		self.ring_pose = ring_pose
		self.ring_obj = ring_obj
		self.protrusion = protrusion
		self.NTF2_bp = Blueprint(data=[ [ x for x in i ] for i in ring_obj.ring_bp.bp_data ] )
		self.NTF2_bp.freeze_all()
		self.NTF2_bp.reindex_blueprint()
		self.generate_bp()
		self.NTF2_dict = {}
		self.has_cHelix = chelix
		self.c_helix_dict = {}
		self.populate_NTF2_dict()

	def write_blueprint(self,fname=None):
		if not fname:
			fname = './NTF2.bp'
		sspairs_str = "SSPAIR "+";".join(self.NTF2_dict['ring_dict']['pairings'])
		self.NTF2_bp.dump_blueprint(fname,[sspairs_str])

	def get_min_fix_res(self):
		'''
		Get the residues to fix during minimization step.
		'''
		R_prev = True
		switch_pos = []
		for n,i in enumerate(self.NTF2_bp.bp_data):
			if i[-1] !='R':
				R_now = False
				if R_now != R_prev:
					switch_pos.append(n+1)
				R_prev = False
			else:
				R_now = True
				if R_now != R_prev:
					switch_pos.append(n)
				R_prev = True
		if self.NTF2_bp.bp_data[-1][-1] != 'R':
			switch_pos.append(len(self.NTF2_bp.bp_data))

		if len(switch_pos)%2 !=0:
			raise ValueError('The intervals for fixing residues are not paired, something went wrong: see: '+','.join(switch_pos))

		return switch_pos

	def get_noPro_posL2(self):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		l3n = bp.segment_dict['L3'].bp_data[1][0]
		l3c = bp.segment_dict['L3'].bp_data[-1][0]
		h3n = bp.segment_dict['H3'].bp_data[0][0]
		return "%d,%d,%d"%(l3n,l3c,h3n)

	def create_helix_csts(self,fname=None):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		cst_string = ''
		cst_string += PerfectHelixCst(bp,1)
		cst_string += PerfectHelixCst(bp,2)
		if not fname:
			fname = './NTF2.csts'
		handle = open(fname,'a')
		handle.write(cst_string)
		handle.close()

	def create_HBond_csts(self,fname=None):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		hbond1_don = bp.segment_dict['L5'].bp_data[0][0]
		hbond1_acc = bp.segment_dict['L3'].bp_data[-1][0]
		cst_string = CircularHBondConstraints(hbond1_don,hbond1_acc)

		hbond2_don = bp.segment_dict['H3'].bp_data[0][0]
		hbond2_acc = bp.segment_dict['L3'].bp_data[-3][0]
		cst_string += CircularHBondConstraints(hbond2_don,hbond2_acc)

		hbond3_don = bp.segment_dict['L3'].bp_data[-1][0]
		hbond3_acc = bp.segment_dict['L3'].bp_data[-4][0]
		cst_string += CircularHBondConstraints(hbond3_don,hbond3_acc)

		hbond4_don = bp.segment_dict['L3'].bp_data[-4][0]
		hbond4_acc_pos = 2*(self.NTF2_dict['ring_dict']["sheet_dict"]["short_arm_l"]) - 1
		hbond4_acc = bp.segment_dict['E6'].bp_data[hbond4_acc_pos][0]
		cst_string += CircularHBondConstraints(hbond4_don,hbond4_acc)

		hbond5_acc = bp.segment_dict['H2'].bp_data[-1][0]
		hbond5_don_pos = 2*(self.NTF2_dict['ring_dict']["sheet_dict"]["short_arm_l"]) - 1
		hbond5_don = bp.segment_dict['E6'].bp_data[hbond4_acc_pos][0]
		cst_string += CircularHBondConstraints(hbond5_don,hbond5_acc)

		hbondH1_acc = bp.segment_dict['H1'].bp_data[-3][0]
		hbondH1_don = bp.segment_dict['L2'].bp_data[0][0]
		cst_string += "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(hbondH1_don,hbondH1_acc)

		if not fname:
			fname = './NTF2.csts'
		handle = open(fname,'a')
		handle.write(cst_string)
		handle.close()

	def create_H1C_csts(self,fname=None):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		E3_res1 = bp.segment_dict['E3'].bp_data[0][0] - 1
		E3_res3 = bp.segment_dict['E3'].bp_data[2][0] - 1
		H1_C = bp.segment_dict['H1'].bp_data[-2][0]
		#H2_key_pos = bp.segment_dict['H2'].bp_data[-7][0]
		H2_key_pos = bp.segment_dict['H2'].bp_data[0][0]
		H3_key_pos = bp.segment_dict['H3'].bp_data[7][0]
		cst_string = ""
		if self.Opening == 'Tropical':
			E3_A_ABEGO = max( [ i[0] for i in bp.segment_dict['E3'].bp_data if i[2] == 'EA' ] )
			E3_res1 = E3_A_ABEGO - 3
			res1_dist = 6
			cst_string = "AtomPair CA %i CA %i HARMONIC %0.2f 1.5\n" %(E3_res1,H1_C,res1_dist)
			# Angle cst to avoid the C term of the helix moving forward too much
			E3_res2 = E3_A_ABEGO - 1
			ang = np.deg2rad(90.0)
			ang_tol=np.deg2rad(20.0)
			cst_string += "Angle CA %i CA %i CA %i CIRCULARHARMONIC %3.1f %3.1f\n" %(E3_res2,E3_res1,H1_C,ang,ang_tol)
			# Also keep H2 N far from H3:
			cst_string += "AtomPair CA %i CA %i HARMONIC 9.0 3.0 \n" %(H2_key_pos,H3_key_pos)
		else:
			if self.NTF2_dict['ring_dict']['sheet_dict']['Second_bulge_E3'] and (self.protrusion <= 6) :
				E3_A_ABEGO = min( [ i[0] for i in bp.segment_dict['E3'].bp_data if i[2] == 'EA' ] )
				E3_res1 = E3_A_ABEGO - 1
				H1_C = bp.segment_dict['H1'].bp_data[-2][0]
				# Knob-hole arrangement between hole on C term of H1 and residue of E3
				cst_string = "AtomPair N %i CB %i HARMONIC 6.0 2.0\n" %(H1_C,E3_res1)
				E3_bulge_in_angle = math.radians(150)
				E3_bulge_in_angle_sd = math.radians(90)
			else:
				E3_res1 = bp.segment_dict['E3'].bp_data[1][0]
				H1_C = bp.segment_dict['H1'].bp_data[-2][0]
				# Knob-hole arrangement between hole on C term of H1 and residue of E3
				cst_string = "AtomPair N %i CB %i HARMONIC 6.0 2.0\n" %(H1_C,E3_res1)
				E3_bulge_in_angle = math.radians(150)
				E3_bulge_in_angle_sd = math.radians(90)
			# Knob-hole arrangement between hole on N term of H2 and residue of H3
			cst_string += "AtomPair CB %i N %i HARMONIC 7.0 3.0\n" %(H2_key_pos,H3_key_pos)
		if not fname:
			fname = './NTF2.csts'
		handle = open(fname,'a')
		handle.write(cst_string)
		handle.close()

	def create_H1N_csts(self,fname=None):

		bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()

		if self.NTF2_dict['ring_dict']['sheet_dict']['short_arm_l'] == 2:
			short_arm_key_pos = bp.segment_dict['E6'].bp_data[1][0]
		else:
			short_arm_key_pos = bp.segment_dict['E6'].bp_data[0][0] - 1

		if self.NTF2_dict['h2_len'] == 7:
			H1_n1 = bp.segment_dict['H1'].bp_data[-11][0]
			H1_n2 = H1_n1 - 3
			#H1_n3 = H1_n1 - 4
		else:
			H1_n1 = bp.segment_dict['H1'].bp_data[-15][0]
			H1_n2 = H1_n1 - 3
			#H1_n3 = H1_n1 - 4

		# Knob-hole arrangement between hole on N term of H1 and residue of loop between E5 and E6.
		cst_string = "AtomPair CA %i N %i HARMONIC 8.0 2.0\n" %(short_arm_key_pos,H1_n2)
		angle = math.radians(130)
		angle_sd = math.radians(90)
		#cst_string += "Angle C %i CA %i N %i HARMONIC %0.2f %0.2f\n" %(short_arm_key_pos,short_arm_key_pos,H1_n2,angle,angle_sd)
		if not fname:
			fname = './NTF2.csts'
		handle = open(fname,'a')
		handle.write(cst_string)
		handle.close()

	def AutomaticCstCreation(self, bp_fname=None, min_fname=None):
		if not bp_fname:
			bp_fname = 'NterHs_bp.csts'
		if not min_fname:
			min_fname = 'NterHs_min.csts'
		self.create_HBond_csts(fname=bp_fname)
		self.create_H1C_csts(fname=bp_fname)
		self.create_H1N_csts(fname=bp_fname)
		##############################
		self.create_helix_csts(fname=min_fname)
		self.create_HBond_csts(fname=min_fname)
		self.create_H1C_csts(fname=min_fname)
		self.create_H1N_csts(fname=min_fname)

	def get_long_arm_inward_pos(self):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		E4 = bp.segment_dict['E4'].bp_data[-1][0]
		E5 = bp.segment_dict['E5'].bp_data[0][0]
		long_arm_length = self.NTF2_dict['ring_dict']['sheet_dict']['long_arm_l']
		positions = [E4,E5]
		for i in range(long_arm_length-1):
			E4 = E4-2
			E5 = E5+2
			positions = positions + [E4,E5]

		return positions
	def get_H1_positions(self):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		H1_all = bp.segment_dict['H1'].bp_data
		positions = [ i[0] for i in H1_all ]
		return positions

	def get_H3_pos(self):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		H3_all = bp.segment_dict['H3'].bp_data
		positions = [ i[0] for i in H3_all ]
		return positions

	def get_H4_positions(self):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		H4_all = bp.segment_dict['H4'].bp_data
		positions = [ i[0] for i in H4_all ]
		return positions

	def get_gly_resc_phe(self):
		NTF2_bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		sheet_dict = self.NTF2_dict['ring_dict']['sheet_dict']
		sheet_obj = CreateSheetObjFromDict(sheet_dict)
		# Here I will return positions in a list, first the ones associated with the main bulge, then the ones with the secondary bulge
		phe_rel_positions = []
		rel_MainE5_phe_positions = sheet_obj.get_MainBulge_hinges()[1]
		phe_rel_positions.append(rel_MainE5_phe_positions)
		if sheet_obj.sheet_data['Second_bulge_E3']:
			rel_SecE5_phe_positions = sheet_obj.get_SecBulge_hinges()[1]
			phe_rel_positions.append(rel_SecE5_phe_positions)
		phe_abs_positions = []
		for i in phe_rel_positions:
			phe_abs_positions.append(NTF2_bp.segment_dict['E5'].bp_data[i][0])
		return phe_abs_positions

	def get_gly_resc_gly(self):
		NTF2_bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		sheet_dict = self.NTF2_dict['ring_dict']['sheet_dict']
		sheet_obj = CreateSheetObjFromDict(sheet_dict)
		# Here I will return positions in a list, first the ones associated with the main bulge, then the ones with the secondary bulge
		gly_rel_positions = []
		rel_MainE4_gly_positions = sheet_obj.get_MainBulge_hinges()[0]
		gly_rel_positions.append(rel_MainE4_gly_positions)
		if sheet_obj.sheet_data['Second_bulge_E3']:
			rel_SecE4_gly_positions = sheet_obj.get_SecBulge_hinges()[0]
			gly_rel_positions.append(rel_SecE4_gly_positions)
		if sheet_obj.sheet_data['CurvedLongArm']:
			gly_rel_positions.append( sheet_obj.get_CurvedLongArm_hinges() )
		gly_abs_positions = []
		for i in gly_rel_positions:
			gly_abs_positions.append(NTF2_bp.segment_dict['E4'].bp_data[i][0])
		return gly_abs_positions
	def get_tight_PRO(self):
		NTF2_bp = Blueprint(data=[ [ x for x in i ] for i in self.NTF2_bp.bp_data ] )
		sheet_dict = self.NTF2_dict['ring_dict']['sheet_dict']
		if sheet_dict["ExtendedE6"]: return []
		sheet_obj = CreateSheetObjFromDict(sheet_dict)
		relative_position = len(NTF2_bp.segment_dict['E5'].bp_data) - len(NTF2_bp.segment_dict['E6'].bp_data)
		absolute_position = NTF2_bp.segment_dict['E5'].bp_data[relative_position][0]
		tight_pro_positions = [ absolute_position-1,absolute_position,absolute_position+1]
		return tight_pro_positions

def CreateBasicNTF2fromDict(NTF2dict,NTF2_is_complete=True,db=None):
	ring_obj = CreateRingObjFromDict(NTF2dict["ring_dict"],db=db)
	loopTwoABEGO = NTF2dict["loopTwoABEGO"]
	Opening = NTF2dict["Opening"]
	h1_len = NTF2dict["h1_len"]
	h2_len = NTF2dict["h2_len"]
	protrusion = NTF2dict["protrusion"]

	if "has_cHelix" in NTF2dict.keys():
		has_cHelix = NTF2dict["has_cHelix"]
	if "chelix" in NTF2dict.keys():
		chelix = NTF2dict["chelix"]

	NTF2_obj = BasicBeNTF2(ring_obj, h1_len = h1_len, loopTwoABEGO = loopTwoABEGO, h2_len = h2_len, Opening=Opening, protrusion=protrusion, chelix=has_cHelix )
	if NTF2_is_complete and has_cHelix:
		# Disposable object just to correctly populate C-term H blueprint
		disposable_obj = SetUpCHelixStep(NTF2_obj)
	return NTF2_obj


def Get_E5Bulge_longArm_dist(ring_pose, ring_obj):

	ring_bp = Blueprint(data=[ [ x for x in i ] for i in ring_obj.ring_bp.bp_data ] )
	ring_bp.reindex_blueprint()

	E6_bulge_pos = [ i[0] for i in ring_bp.segment_dict['E6'].bp_data if i[2] == 'EA' ][0]
	E3_N_pos = ring_bp.segment_dict['E3'].bp_data[1][0]
	H3_C_pos = ring_bp.segment_dict['H2'].bp_data[-3][0]

	E6_bulge_CA = ring_pose.residue(E6_bulge_pos).xyz('CA')
	E3_N_CA = ring_pose.residue(E3_N_pos).xyz('CA')
	H3_C_CA = ring_pose.residue(H3_C_pos).xyz('CA')
	#delta = ( E3_N_CA - H3_C_CA )
	#midpoint = numeric.xyzVector_double_t(delta[0]/2, delta[1]/2, delta[2]/2 )
	midpoint = numeric.xyzVector_double_t( (E3_N_CA[0]+H3_C_CA[0])/2, (E3_N_CA[1]+H3_C_CA[1])/2, (E3_N_CA[2]+H3_C_CA[2])/2 )
	distance_v = midpoint - E6_bulge_CA
	distance = distance_v.length()
	return distance

def Get_CrossSheet_dist(ring_pose, ring_obj):
	ring_bp = Blueprint(data=[ [ x for x in i ] for i in ring_obj.ring_bp.bp_data ] )
	ring_bp.reindex_blueprint()

	#E3_N_pos = ring_bp.segment_dict['E3'].bp_data[0][0]
	H3_C_pos = ring_bp.segment_dict['H2'].bp_data[-1][0]
	E3_C_pos = ring_bp.segment_dict['E3'].bp_data[-1][0]

	#E3_N_CA = ring_pose.residue(E3_N_pos).xyz('CA')
	H3_C_CA = ring_pose.residue(H3_C_pos).xyz('CA')
	E3_C_CA = ring_pose.residue(E3_C_pos).xyz('CA')

	#distance_v = E3_N_CA - E3_C_CA
	distance_v = H3_C_CA - E3_C_CA
	distance = distance_v.length()
	return distance

def Get_HP_longArm_dist(ring_pose, ring_obj):
	ring_bp = Blueprint(data=[ [ x for x in i ] for i in ring_obj.ring_bp.bp_data ] )
	ring_bp.reindex_blueprint()

	#E3_N_pos = ring_bp.segment_dict['E3'].bp_data[0][0]
	E4_C_pos = ring_bp.segment_dict['E4'].bp_data[-1][0]
	E1_C_pos = ring_bp.segment_dict['E1'].bp_data[-1][0]

	#E3_N_CA = ring_pose.residue(E3_N_pos).xyz('CA')
	E4_C_CA = ring_pose.residue(E4_C_pos).xyz('CA')
	E1_C_CA = ring_pose.residue(E1_C_pos).xyz('CA')

	#distance_v = E3_N_CA - E3_C_CA
	distance_v = E1_C_CA - E4_C_CA
	distance = distance_v.length()
	return distance

def Get_E6_C_E5_N_dist(ring_pose, ring_obj):
	ring_bp = Blueprint(data=[ [ x for x in i ] for i in ring_obj.ring_bp.bp_data ] )
	ring_bp.reindex_blueprint()

	#E3_N_pos = ring_bp.segment_dict['E3'].bp_data[0][0]
	E5_N_pos = ring_bp.segment_dict['E5'].bp_data[0][0]
	E6_C_pos = ring_bp.segment_dict['E6'].bp_data[-1][0]

	#E3_N_CA = ring_pose.residue(E3_N_pos).xyz('CA')
	E5_N_CA = ring_pose.residue(E5_N_pos).xyz('CA')
	E6_C_CA = ring_pose.residue(E6_C_pos).xyz('CA')

	#distance_v = E3_N_CA - E3_C_CA
	distance_v = E5_N_CA - E6_C_CA
	distance = distance_v.length()
	return distance

def GetShortArm_H3_7_dist(ring_pose, ring_obj):
	ring_bp = Blueprint(data=[ [ x for x in i ] for i in ring_obj.ring_bp.bp_data ] )
	ring_bp.reindex_blueprint()

	#E3_N_pos = ring_bp.segment_dict['E3'].bp_data[0][0]
	H3_key_pos = ring_bp.segment_dict['H2'].bp_data[7][0]
	if ring_obj.ring_dict['sheet_data']['short_arm_l'] == 2:
		short_arm_key_pos = ring_bp.segment_dict['E6'].bp_data[1][0]
	else:
		short_arm_key_pos = ring_bp.segment_dict['E6'].bp_data[0][0] - 1

	H3_key_pos_CA = ring_pose.residue(H3_key_pos).xyz('CA')
	short_arm_key_pos_CA = ring_pose.residue(short_arm_key_pos).xyz('CA')

	distance_v = H3_key_pos_CA - short_arm_key_pos_CA
	distance = distance_v.length()
	return distance

def TP_is_viable(TP_pose, BeNTF2_obj):
	BeNTF2_bp = Blueprint(data=[ [ x for x in i ] for i in BeNTF2_obj.NTF2_bp.bp_data ] )
	BeNTF2_bp.reindex_blueprint()

	#E3_N_pos = ring_bp.segment_dict['E3'].bp_data[0][0]
	H2_N_pos = BeNTF2_bp.segment_dict['H2'].bp_data[0][0]
	H3_C_pos = BeNTF2_bp.segment_dict['H3'].bp_data[-1][0]

	#E3_N_CA = ring_pose.residue(E3_N_pos).xyz('CA')
	H2_N_CA = TP_pose.residue(H2_N_pos).xyz('CA')
	H3_C_CA = TP_pose.residue(H3_C_pos).xyz('CA')

	#distance_v = E3_N_CA - E3_C_CA
	distance_v = H2_N_CA - H3_C_CA
	distance = distance_v.length()
	if distance < 12.0:
		return False
	else:
		return True

def RingCanBeTP_NTF2(ring_pose,ring_obj):
	mouth_size = Get_HP_longArm_dist(ring_pose, ring_obj)
	dist = Get_E5Bulge_longArm_dist(ring_pose, ring_obj)
	if ( mouth_size < 13.0 ) or ( mouth_size > 16.0 ): #or ( dist < 22 ) :
			print("This ring has either a mouth that is too wide or too small, reverting to Classic opening with mouth_size = %0.2f"%mouth_size)
			return False
	else:
			return True

def SetUpNtermStep(ring_pose, ring_obj, tropical = False, chelix = False):
	protrusion = Get_RinglongArm_protrusionDist(ring_pose,ring_obj)
	H3_E6_B_dist = Get_E5Bulge_longArm_dist(ring_pose, ring_obj)
	E6_C_E5_N = Get_E6_C_E5_N_dist(ring_pose, ring_obj)
	print('H3_E6_B_dist: %0.4f'%H3_E6_B_dist)
	cross_sheet_dist = Get_CrossSheet_dist(ring_pose, ring_obj)
	if tropical:
		#TestTropical(ring_pose, ring_obj)
		Opening_type = 'Tropical'
		mouth_size = Get_HP_longArm_dist(ring_pose, ring_obj)
		if ( mouth_size < 13.0 ) or ( mouth_size > 16.0 ):
			print("This ring has either a mouth that is too wide or too small, reverting to Classic opening")
			tropical = False
		else:
			chelix = True
			H2_len = 7
			H1_length = 14

	if not tropical:
		Opening_type = 'Classic'
		# C helix check:
		if (E6_C_E5_N > 18.5) or (E6_C_E5_N < 15.0):
			print("E6_C_E5_N: %0.2f"%E6_C_E5_N)
			chelix = False
			print("The mouth is too small or too big for a C-term helix")

		if (H3_E6_B_dist > 25):
			H2_len = 11
		else:
			H2_len = 7
		print('cross_sheet_dist: %0.4f'%cross_sheet_dist)
		#H1_length = int( round( (cross_sheet_dist/5.4)*3.6 - 3.0) )

		H1_length = 16 + ( H2_len-7 ) + 3

	return BasicBeNTF2(ring_obj, ring_pose = ring_pose , h1_len = H1_length, loopTwoABEGO = 'GB', h2_len = H2_len, Opening=Opening_type, protrusion=protrusion , chelix = chelix)

class NTF2_CTermHelix():

	def build_bp(self,fname = None):
		if not self.is_tropical_pitcher:
			bp = self.basic_NTF2_obj.NTF2_bp
			bp.reindex_blueprint()
			bp.freeze_all()
			bp.bp_data[-2][-1] = 'R'
			bp.bp_data[-1][-1] = 'R'
			bp.bp_data[-1][-2] = 'L%s'%(self.loopABEGO)
			for i in range(self.h_len):
				bp.bp_data.append([0,'A','HA','R'])
			bp.bp_data.append([0,'A','LX','R'])
		else:
			bp = self.basic_NTF2_obj.NTF2_bp
			bp.reindex_blueprint()
			bp.freeze_all()
			bp.bp_data[-2][-1] = 'R'
			bp.bp_data[-1][-1] = 'R'
			bp.bp_data[-1][-2] = 'LA'
			bp.bp_data.append([0,'A','LB','R'])
			for i in range(self.h_len):
				bp.bp_data.append([0,'A','HA','R'])
			bp.bp_data.append([0,'A','LX','R'])
		return self.basic_NTF2_obj.NTF2_bp

	def write_bp(self,fname = None):
		if not fname:
			fname = './CH_NTF2.bp'
		sspairs_str = "SSPAIR "+";".join(self.basic_NTF2_obj.NTF2_dict['ring_dict']['pairings'])
		self.basic_NTF2_obj.NTF2_bp.dump_blueprint(fname,[sspairs_str])

	def get_no_pro(self):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.basic_NTF2_obj.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		pos = bp.segment_dict['L10'].bp_data[-1][0]
		return pos

	def create_csts(self,fname = None):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.basic_NTF2_obj.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		cst_string = ''
		if not self.is_tropical_pitcher:
			HC_C_term = bp.bp_data[-2][0]
			HC_preC_term = bp.bp_data[-3][0]
			E4_C = bp.segment_dict['E4'].bp_data[-1][0]
			E5_N = E4_C + 3
			cst_string += "AtomPair CA %i CA %i HARMONIC 7.0 3.0\n" %(HC_C_term,E4_C)
			cst_string += "AtomPair CA %i CA %i HARMONIC 7.0 3.0\n" %(HC_preC_term,E4_C)
			cst_string += "AtomPair CA %i CA %i HARMONIC 7.0 3.0\n" %(HC_C_term,E5_N)
			cst_string += "AtomPair CA %i CA %i HARMONIC 7.0 3.0\n" %(HC_preC_term,E5_N)
		else:
			HC_C_term = bp.segment_dict['H4'].bp_data[6][0]
			E2_N = bp.segment_dict['E2'].bp_data[-4][0]
			E1_C = bp.segment_dict['E1'].bp_data[3][0]
			E4_C = bp.segment_dict['E4'].bp_data[-1][0]
			#CSTs for the middle of C-helix
			E1_C_CA = self.basic_NTF2_pose.residue(E1_C).xyz('CA')
			E4_C_CA = self.basic_NTF2_pose.residue(E4_C).xyz('CA')
			midpoint = numeric.xyzVector_double_t( (E1_C_CA[0]+E4_C_CA[0])/2, (E1_C_CA[1]+E4_C_CA[1])/2, (E1_C_CA[2]+E4_C_CA[2])/2 )
			distance_v = midpoint - E1_C_CA
			distance = distance_v.length() - 1.0
			cst_string += "AtomPair CA %i CA %i HARMONIC %0.2f 2.0\n" %(HC_C_term,E2_N,distance)
			cst_string += "AtomPair CA %i CA %i HARMONIC %0.2f 2.0\n" %(HC_C_term,E4_C,distance)
			#CSTs for the C term:
			H3_middle = bp.segment_dict['H3'].bp_data[6][0]
			HC_C_end_a = bp.segment_dict['H4'].bp_data[-1][0]
			distance = 9.0
			dev = 3.0
			HC_C_end_b = HC_C_end_a - 1
			HC_C_end_c = HC_C_end_b - 1
			cst_string += "AtomPair CA %i CA %i HARMONIC %0.2f %0.2f\n" %(HC_C_end_a,H3_middle,distance,dev)
			cst_string += "AtomPair CA %i CA %i HARMONIC %0.2f %0.2f\n" %(HC_C_end_b,H3_middle,distance,dev)
			cst_string += "AtomPair CA %i CA %i HARMONIC %0.2f %0.2f\n" %(HC_C_end_c,H3_middle,distance,dev)
			# hbond csts
			hbond1_don = bp.segment_dict['L10'].bp_data[-1][0]
			pairings = self.basic_NTF2_obj.NTF2_dict['ring_dict']['sheet_dict']['pairings']
			pairing_5_6 = pairings["SSPAIR_5_6"]
			hbond1_acc_relative = pairings["SSPAIR_5_6"] - 1
			hbond1_acc = bp.segment_dict['E5'].bp_data[hbond1_acc_relative][0]
			cst_string += CircularHBondConstraints(hbond1_don,hbond1_acc)
		# Straight helix csts
		cst_string += PerfectHelixCst(bp,4)
		#return cst_string
		self.cst_lines = cst_string

	def write_csts(self,fname = None):
		self.create_csts()
		if not fname:
			fname = './CH_NTF2.csts'
		handle = open(fname,'a')
		handle.write(self.cst_lines)
		handle.close()

	def populate_Ch_dict(self):
		ch_dict = {
				'h1_len':self.h_len, \
				'loopABEGO':self.loopABEGO,\
				'is_tropical_pitcher':self.is_tropical_pitcher,\
				'cst_lines':self.cst_lines \
				}
		return ch_dict

	def __init__(self, basic_NTF2_obj, basic_NTF2_pose = None , h_len = 8, loopABEGO = 'B' ):
		self.loopABEGO = loopABEGO
		self.h_len = h_len
		self.basic_NTF2_pose = basic_NTF2_pose
		self.basic_NTF2_obj = basic_NTF2_obj
		self.basic_NTF2_obj.NTF2_dict['has_cHelix'] = True
		self.is_tropical_pitcher = False
		if self.basic_NTF2_obj.NTF2_dict['Opening'] == 'Tropical':
			self.is_tropical_pitcher = True
		self.bp_w_ch = self.build_bp()
		self.cst_lines = ''
		#self.create_csts()
		self.ch_dict = self.populate_Ch_dict()
		self.basic_NTF2_obj.NTF2_dict['c_helix_dict'] = self.ch_dict

	def get_min_fix_res(self):
		bp = Blueprint(data=[ [ x for x in i ] for i in self.basic_NTF2_obj.NTF2_bp.bp_data ] )
		bp.reindex_blueprint()
		switch_pos = bp.segment_dict['E6'].bp_data[-3][0]
		return switch_pos

def NTF2CanHaveCTermH(BasicNTF2_pose, BasicNTF2_obj=None,ring_obj=None):
	ring_bp = None
	if BasicNTF2_obj:
		ring_bp = Blueprint(data=[ [ x for x in i ] for i in BasicNTF2_obj.NTF2_bp.bp_data ] )
	elif ring_obj:
		ring_bp = Blueprint(data=[ [ x for x in i ] for i in ring_obj.ring_bp.bp_data ] )
	ring_bp.reindex_blueprint()

	#E3_N_pos = ring_bp.segment_dict['E3'].bp_data[0][0]
	E5_N_pos = ring_bp.segment_dict['E5'].bp_data[0][0]
	E6_C_pos = ring_bp.segment_dict['E6'].bp_data[-1][0]

	#E3_N_CA = ring_pose.residue(E3_N_pos).xyz('CA')
	E5_N_CA = BasicNTF2_pose.residue(E5_N_pos).xyz('CA')
	E6_C_CA = BasicNTF2_pose.residue(E6_C_pos).xyz('CA')

	#distance_v = E3_N_CA - E3_C_CA
	distance_v = E5_N_CA - E6_C_CA
	E6_C_E5_N = distance_v.length()

	if (E6_C_E5_N > 18.5) or (E6_C_E5_N < 15.0):
		print("E6_C_E5_N: %0.2f"%E6_C_E5_N)
		print("The mouth is too small or too big for a C-term helix with E6_C_E5_N = %0.2f"%E6_C_E5_N)
		return False
	else:
		return True

def SetUpCHelixStep(BasicNTF2_obj, BasicNTF2_pose = None ):
	if BasicNTF2_obj.NTF2_dict["Opening"] == "Tropical":
		return NTF2_CTermHelix(BasicNTF2_obj, basic_NTF2_pose = BasicNTF2_pose, h_len = 11, loopABEGO = 'AB' )
	else:
		return NTF2_CTermHelix(BasicNTF2_obj, basic_NTF2_pose = BasicNTF2_pose, h_len = 8 , loopABEGO = 'B')


if __name__ == "__main__":
	print("This is the NTF2 toolkit by Benjamin Basanta - June 2017")
