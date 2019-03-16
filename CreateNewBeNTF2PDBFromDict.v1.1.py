#!/software/conda/envs/pyrosetta/bin/python3.7

from CreateBeNTF2_backbone import *
#from BeNTF2_toolkit import *
#import json
#import argparse
#import copy
#import os

'''
def CreateSheetPoseFromObj(sheet_obj,POSE,n_trials_sheet=25,sheet_type=None):
        failed_mover_statuses_names = [ "FAIL_RETRY", "FAIL_DO_NOT_RETRY", "FAIL_BAD_INPUT"]
        failed_mover_statuses = [ protocols.moves.mstype_from_name( name ) for name in failed_mover_statuses_names ]

        init_pose_len = POSE.size()
        print('Attempting to create sheet')
        shee_angle_fname = './%ssheet_angles_%d.csts'%(prefix,trial)
        sheet_bp_fname = './%ssheet_bp_%d.bp'%(prefix,trial)
        sheet = sheet_obj
        sheet.write_cstfile(shee_angle_fname)
        sspairs_str = "SSPAIR "+";".join(sheet.pairings)
        bp = Blueprint( data=[ [ x for x in i ] for i in sheet.blueprint.bp_data ] )
        bp.remodel_all()
        bp.index_to_zero_all()
        bp.bp_data[0][0] = 1
        bp.bp_data[1][0] = 2
        bp.dump_blueprint(sheet_bp_fname,[sspairs_str])

        replaces = ['remodel_cen=%s/remodel_cen.wts'%(db),'fldsgn_cen=%s/fldsgn_cen.wts'%(db),\
                        'cst_file=%s'%shee_angle_fname,'bp=%s'%sheet_bp_fname,'sheet_attempts=%d'%n_trials_sheet]
        tomponents_dict = sheet.get_Tomponents_flags_dict()
        print( tomponents_dict )
        for key in tomponents_dict.keys():
                replaces.append("%s=%s"%(key,tomponents_dict[key]))
        print('Creating sheet mover')
        mover_from_XML = MakeMoverFromXML('%s/CreateSheet.xml'%db,replaces,POSE)
        print('Applying')
        try:mover_from_XML.apply(POSE)
        except RuntimeError:
                print('Failed to create create sheet')
                os.remove('./%s'%shee_angle_fname)
                os.remove('./%s'%sheet_bp_fname)
                return False

        else:
                if mover_from_XML.get_last_move_status() in failed_mover_statuses :
                        print('Failed to create create sheet')
                        os.remove('./%s'%shee_angle_fname)
                        os.remove('./%s'%sheet_bp_fname)
                        return False
                else:
                        os.remove('./%s'%shee_angle_fname)
                        os.remove('./%s'%sheet_bp_fname)
                        if init_pose_len == POSE.size():
                                print('Failed to create create sheet')
                                return False
                        else:
                                print('Successfully created sheet')
                                ### EXTRACT FEATURES: ###
                                dist = Get_SheetE6bulge_E3N_dist(POSE,sheet_obj)
                                long_arm_ang = Get_SheetBase_longArm_angle(POSE,sheet_obj)
                                shoot = Get_longArm_shootingAngle(POSE,sheet_obj)
                                prot = Get_longArm_protrusionDist(POSE,sheet_obj)
                                sheet_type = GetSheetType(sheet_obj,db)
                                base = sheet_obj.sheet_data['base_width']
                                long_arm = sheet_obj.sheet_data['long_arm_l']
                                short_arm = sheet_obj.sheet_data['short_arm_l']
                                sec_bulge = 0 if not sheet_obj.sheet_data["Second_bulge_E3"] else 1
                                main_bulge_curve = sheet_obj.sheet_data["E3_MainBulgeCurve"]
                                features_v = { 'dist':dist,'long_arm_ang':long_arm_ang,'shoot':shoot,\
                                               'prot':prot,'sheet_type':sheet_type,'base':base,'long_arm':long_arm, \
                                               'short_arm':short_arm,'sec_bulge':sec_bulge, \
                                               'main_bulge_curve':main_bulge_curve }

                                dict_string = json.dumps(sheet_obj.sheet_data)
                                add_BeNTF2_dict_comment(POSE,dict_string)
                                core.pose.add_comment(POSE,"FEATURES",json.dumps(features_v))
                                return True

def AssembleRing(POSE,Ring,ring_trials=10):
        failed_mover_statuses_names = [ "FAIL_RETRY", "FAIL_DO_NOT_RETRY", "FAIL_BAD_INPUT"]
        failed_mover_statuses = [ protocols.moves.mstype_from_name( name ) for name in failed_mover_statuses_names ]
        sheet_len = POSE.size()
        ring_bp_name = './%sring_%d.bp'%(prefix,trial)
        ring_allCST_name = './%sring_CST_%d.cst'%(prefix,trial)
        Ring.write_blueprint(ring_bp_name)
        Ring.dump_all_csts(ring_allCST_name)
        replaces = ['remodel_cen=%s/remodel_cen.wts'%(db),'all_cst_file=%s'%ring_allCST_name,'bp=%s'%ring_bp_name,'ring_attempts=%d'%ring_trials]
        fix_res = Ring.get_min_fix_res()
        frontal_HP_res = Ring.FrontHPHbonds()
        replaces.append("fix_min=%d-%d"%(fix_res[0],fix_res[1]))
        replaces.append("no_pro=%s"%Ring.get_noPro_posL2())
        replaces.append("res1=%d"%frontal_HP_res[0])
        replaces.append("res2=%d"%frontal_HP_res[1])
        mover_from_XML = MakeMoverFromXML('%s/RingConnection.xml'%db,replaces,POSE)
        try:mover_from_XML.apply(POSE)
        except RuntimeError:
                        print('Unable to create this particular ring')
                        os.remove('./%s'%ring_bp_name)
                        os.remove('./%s'%ring_allCST_name)
                        #os.remove('./%s'ring_dict_name)
                        return False
        else:
                if sheet_len == POSE.size() or (mover_from_XML.get_last_move_status() in failed_mover_statuses ):
                        print('Unable to create this particular ring')
                        os.remove('./%s'%ring_bp_name)
                        os.remove('./%s'%ring_allCST_name)
                        return False
                else:
                        os.remove('./%s'%ring_bp_name)
                        os.remove('./%s'%ring_allCST_name)
                        return True

def AssembleBasicBeNTF2(POSE,BasicBeNTF2_obj,trials=10):
        failed_mover_statuses_names = [ "FAIL_RETRY", "FAIL_DO_NOT_RETRY", "FAIL_BAD_INPUT"]
        failed_mover_statuses = [ protocols.moves.mstype_from_name( name ) for name in failed_mover_statuses_names ]
        dict_string = json.dumps(BasicBeNTF2_obj.NTF2_dict)
        ring_len = POSE.size()
        BeNTF2_bp_name = './%sBasicBeNTF2_%d.bp'%(prefix,trial)
        BeNTF2_bpCST_name = './%sBasicBeNTF2_bpCST_%d.cst'%(prefix,trial)
        BeNTF2_minCST_name = './%sBasicBeNTF2_minCST_%d.cst'%(prefix,trial)
        #ring_hbCST_name = './%sring_hbCST_%d.bp'%(prefix,trial)
        BasicBeNTF2_obj.write_blueprint(BeNTF2_bp_name)
        BasicBeNTF2_obj.AutomaticCstCreation(bp_fname=BeNTF2_bpCST_name,min_fname=BeNTF2_minCST_name)
        #Ring.dump_Hbonds_csts(ring_hbCST_name)
        replaces = ['remodel_cen=%s/remodel_cen.wts'%(db),'all_cst_file=%s'%BeNTF2_minCST_name,'bp=%s'%BeNTF2_bp_name,'BasicBeNTF2_attempts=%d'%trials]
        fix_res = BasicBeNTF2_obj.get_min_fix_res()
        #frontal_HP_res = Ring.FrontHPHbonds()
        replaces.append("fix_min=%d-%d"%(fix_res[0],fix_res[1]))
        replaces.append("no_pro=%s"%BasicBeNTF2_obj.get_noPro_posL2())
        #replaces.append("res1=%d"%frontal_HP_res[0])
        #replaces.append("res2=%d"%frontal_HP_res[1])
        mover_from_XML = MakeMoverFromXML('%s/basic_BeNTF2.xml'%db,replaces,POSE)
        try:mover_from_XML.apply(POSE)
        except RuntimeError:
                print('Unable to create this N term helices')
                os.remove(BeNTF2_bp_name)
                os.remove(BeNTF2_bpCST_name)
                os.remove(BeNTF2_minCST_name)
                return (dict_string,False)
        else:
                if ring_len == POSE.size() or (mover_from_XML.get_last_move_status() in failed_mover_statuses ) :
                        print('Unable to create N term helices')
                        os.remove(BeNTF2_bp_name)
                        os.remove(BeNTF2_bpCST_name)
                        os.remove(BeNTF2_minCST_name)
                        return False
                else:
                        os.remove(BeNTF2_bp_name)
                        os.remove(BeNTF2_bpCST_name)
                        os.remove(BeNTF2_minCST_name)
                        return True

def AddCtermHelix(BasicBeNTF2_obj,BasicNTF2_pose=None,trials=25):
        failed_mover_statuses_names = [ "FAIL_RETRY", "FAIL_DO_NOT_RETRY", "FAIL_BAD_INPUT"]
        failed_mover_statuses = [ protocols.moves.mstype_from_name( name ) for name in failed_mover_statuses_names ]
        basic_len = BasicNTF2_pose.size()
        CH_BeNTF2_bp_name = './%sCH_BasicBeNTF2_%d.bp'%(prefix,trial)
        CH_BeNTF2_bpCST_name = './%sCH_BasicBeNTF2_bpCST_%d.cst'%(prefix,trial)
        CH_BeNTF2_obj = BeNTF2_obj['c_helix_dict']
        CH_BeNTF2_obj.write_csts(fname = CH_BeNTF2_bpCST_name)
        CH_BeNTF2_obj.write_bp(fname = CH_BeNTF2_bp_name)
        no_pro = CH_BeNTF2_obj.get_no_pro()
        replaces = ['remodel_cen=%s/remodel_cen.wts'%(db),'all_cst_file=%s'%CH_BeNTF2_bpCST_name,'bp=%s'%CH_BeNTF2_bp_name,'CH_BeNTF2_attempts=%d'%trials,'no_pro=%d'%no_pro]
        fix_res = CH_BeNTF2_obj.get_min_fix_res()
        replaces.append("fix_min=1-%d"%(fix_res))
        mover_from_XML = MakeMoverFromXML('%s/CH_BeNTF2.xml'%db,replaces,BasicNTF2_pose)
        try:mover_from_XML.apply(BasicNTF2_pose)
        except RuntimeError:
                print('Unable to create this C term helix')
                os.remove(CH_BeNTF2_bp_name)
                os.remove(CH_BeNTF2_bpCST_name)
                return False
        else:
                if basic_len == POSE.size() or (mover_from_XML.get_last_move_status() in failed_mover_statuses ) :
                        print('Unable to create this C term helix')
                        os.remove(CH_BeNTF2_bp_name)
                        os.remove(CH_BeNTF2_bpCST_name)
                        return False
                else:
                        os.remove(CH_BeNTF2_bp_name)
                        os.remove(CH_BeNTF2_bpCST_name)
                        return True
'''
def CreateBeNTF2PoseFromDict(BeNTF2dict,POSE,arch_len=None):
        init_pose_len = POSE.size()
        ring_data_dict = BeNTF2dict["ring_dict"]
        sheet_dict = BeNTF2dict["ring_dict"]["sheet_dict"]
        #sheet_obj = CreateSheetObjFromDict(sheet_dict)
        #success = CreateSheetPoseFromObj(sheet_obj,POSE,n_trials_sheet=25)

        sheet_obj,success = CreateSheet(POSE,n_trials_sheet=25,sheet_type=ring_data_dict["sheet_type"],\
                              prefix=prefix,trial=trial,db=db,options=options,arch_len=arch_len)
        if not success: return sheet_obj,False
        Base_LArm_angle =  Get_SheetBase_longArm_angle(POSE,sheet_obj)
        Dist_to_bulgeE6 = Get_SheetE6bulge_E3N_dist(POSE,sheet_obj)
        if sheet_obj.sheet_data['CurvedLongArm']:
                if sheet_obj.sheet_data["long_arm_l"] == 3:
                        if Base_LArm_angle > 60:
                                print('Sheet is outward-facing, long arm == 3 and Base_LArm_angle > 60, discarding')
                                return sheet_obj,False
                                #raise Exception('Sheet is outward-facing, discarding')
                if sheet_obj.sheet_data["long_arm_l"] == 4:
                        if Base_LArm_angle > 50:
                                print('Sheet is outward-facing, long arm == 4 and Base_LArm_angle > 50, discarding')
                                return sheet_obj,False
                                #raise Exception('Sheet is outward-facing, discarding')
                        elif Dist_to_bulgeE6 > 29:
                                print('Sheet is outward-facing, long arm == 4 and Dist_to_bulgeE6 > 29, discarding')
                                return sheet_obj,False

        Ring = CreateRingObjFromDict(ring_data_dict,db=db)
        success = AssembleRing(POSE,Ring,ring_trials=10,\
                              prefix=prefix,trial=trial,db=db,options=options)

        if not success: return Ring,False

        can_have_cHelix = NTF2CanHaveCTermH(POSE, ring_obj=Ring)

        if BeNTF2_obj.has_cHelix and not can_have_cHelix:
            success = False

        BeNTF2_obj = CreateBasicNTF2fromDict(BeNTF2dict,NTF2_is_complete=False,db=db)
        success = AssembleBasicBeNTF2(POSE,BeNTF2_obj,trials=25,\
                                      prefix=prefix,trial=trial,db=db,options=options)

        if not success: return BeNTF2_obj,False

        if BeNTF2_obj.has_cHelix:
                success = AddCtermHelix(BeNTF2_obj,BasicNTF2_pose=POSE,trials=25,\
                                        prefix=prefix,trial=trial,db=db,options=options)

        if success: return BeNTF2_obj,True
        else: return BeNTF2_obj,False
'''
def add_SS_positions(POSE,NTF2_obj):
        info = POSE.pdb_info()
        H1_label = 'H1pos'
        H1_positions = NTF2_obj.get_H1_positions()
        for i in H1_positions:
                info.add_reslabel(i,H1_label)

        if NTF2_obj.NTF2_dict["has_cHelix"]:
                H4_label = 'H4pos'
                H4_positions = NTF2_obj.get_H4_positions()
                for i in H4_positions:
                        info.add_reslabel(i,H4_label)

        long_HP_label = 'longHPpos'
        long_HP_positions = NTF2_obj.get_long_arm_inward_pos()
        for i in long_HP_positions:
                info.add_reslabel(i,long_HP_label)
'''
if __name__ == "__main__":
        ########## Option system: ###########
        argparser = argparse.ArgumentParser(description='Create a BeNTF2 backbone.')
        argparser.add_argument('-database', type=str,help='XML file to use for ring generation, second step')
        argparser.add_argument('-input_pdb', type=str,help='sheet pdb')
        argparser.add_argument('-nstruct', type=int,help='Max number of structures to output. May output less if some of them fail.')
        argparser.add_argument('-prefix', type=str,help='Prefix to add to output names.')
        args = argparser.parse_args()
        prefix=''
        if args.prefix:
                prefix = args.prefix
        else:
                prefix = os.path.basename(args.input_pdb).split('_')[0]
        print("Using prefix: %s"%prefix)

        db = args.database
        general_flags = ['-rama_map %s/Rama_XPG_3level.txt'%db,\
                '-picking_old_max_score 1',\
                '-restore_talaris_behavior',
                '-out:file:pdb_comments']# ,\
                #'-mute all']

        init( extra_options=" ".join(general_flags) )
        options = basic.options.process()

        PDB_fname = args.input_pdb
        pdb_handle = open(PDB_fname,'r')
        lines = [ line[:-1] for line in pdb_handle.readlines()]
        dict_line = [ line for line in lines if 'BENTF2DICT' in line ][0]
        features_lines = [ line for line in lines if 'FEATURES' in line ]
        pdb_handle.close()
        dict_string = ' '.join(dict_line.split()[1:])
        BeNTF2dict = json.loads(dict_string)
        features_v = []
        for trial in range(args.nstruct):
                if os.path.isfile('%s_%d.pdb'%(prefix,trial)):
                        continue
                success =  False
                for attempt in range(20):
                        POSE = pose_from_file( filename='%s/input.pdb'%db )
                        if len(features_lines)==1:
                            feat_string = ' '.join(features_lines[0].split()[1:])
                            features_dict = json.loads(feat_string)
                            arch_len_val = features_dict["dist"]
                            BeNTF2_obj,success = CreateBeNTF2PoseFromDict(BeNTF2dict,POSE,arch_len=arch_len_val)
                        else:
                            BeNTF2_obj,success = CreateBeNTF2PoseFromDict(BeNTF2dict,POSE)
                        if success: break
                if not success:
                        print("Unable to make this NTF2 again in the number of iterations provided... sorry...")
                        continue

                # I decided not to add TRP lock and gly positions, those can be automatically detected later:

                POSE.data().clear()
                sheet_dict = BeNTF2dict["ring_dict"]["sheet_dict"]
                sheet_obj = NTF2_sheet(**sheet_dict)
                dist = Get_SheetE6bulge_E3N_dist(POSE,sheet_obj)
                long_arm_ang = Get_SheetBase_longArm_angle(POSE,sheet_obj)
                shoot = Get_longArm_shootingAngle(POSE,sheet_obj)
                prot = Get_longArm_protrusionDist(POSE,sheet_obj)
                sheet_type = GetSheetType(sheet_obj,db)
                base = sheet_obj.sheet_data['base_width']
                long_arm = sheet_obj.sheet_data['long_arm_l']
                short_arm = sheet_obj.sheet_data['short_arm_l']
                sec_bulge = 0 if not sheet_obj.sheet_data["Second_bulge_E3"] else 1
                main_bulge_curve = sheet_obj.sheet_data["E3_MainBulgeCurve"]
                features_v = { 'dist':dist,'long_arm_ang':long_arm_ang,'shoot':shoot,\
                				'prot':prot,'sheet_type':sheet_type,'base':base,'long_arm':long_arm,\
                				'short_arm':short_arm,'sec_bulge':sec_bulge,'main_bulge_curve':main_bulge_curve }

                dict_string = json.dumps(BeNTF2_obj.NTF2_dict)
                add_SS_labels(POSE,BeNTF2_obj)
                add_BeNTF2_dict_comment(POSE,dict_string)
                pocket_positions = detect_pocket_residues(POSE,BeNTF2_obj)
                add_pocket_labels(POSE,pocket_positions)
                add_SS_positions(POSE,BeNTF2_obj)
                core.pose.add_comment(POSE,"FEATURES",json.dumps(features_v))

                try: POSE.dump_file( '%s_%d.pdb'%(prefix,trial) )
                except RuntimeError:
                        print('Unable to create file, skipping')
                print("Successfully recreated BeNTF2!")
