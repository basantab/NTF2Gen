#!/software/conda/envs/pyrosetta/bin/python3.7

from CreateBeNTF2_backbone import *

def CreateBeNTF2PoseFromDict(BeNTF2dict,POSE,arch_len=None):
        init_pose_len = POSE.size()
        ring_data_dict = BeNTF2dict["ring_dict"]
        sheet_dict = BeNTF2dict["ring_dict"]["sheet_dict"]
        features_v = {}
        print("Beginning Sheet construction.")
        sys.stdout.flush()
        sheet_obj,success = CreateSheet(POSE,n_trials_sheet=25,sheet_type=ring_data_dict["sheet_type"],\
                              prefix=prefix,trial=trial,db=db,options=options,arch_len=arch_len)
        if not success: return features_v,sheet_obj,False
        Base_LArm_angle =  Get_SheetBase_longArm_angle(POSE,sheet_obj)
        Dist_to_bulgeE6 = Get_SheetE6bulge_E3N_dist(POSE,sheet_obj)
        if sheet_obj.sheet_data['CurvedLongArm']:
                if sheet_obj.sheet_data["long_arm_l"] == 3:
                        if Base_LArm_angle > 60:
                                print('Sheet is outward-facing, long arm == 3 and Base_LArm_angle > 60, discarding')
                                return features_v,sheet_obj,False
                                #raise Exception('Sheet is outward-facing, discarding')
                if sheet_obj.sheet_data["long_arm_l"] == 4:
                        if Base_LArm_angle > 50:
                                print('Sheet is outward-facing, long arm == 4 and Base_LArm_angle > 50, discarding')
                                return features_v,sheet_obj,False
                                #raise Exception('Sheet is outward-facing, discarding')
                        elif Dist_to_bulgeE6 > 29:
                                print('Sheet is outward-facing, long arm == 4 and Dist_to_bulgeE6 > 29, discarding')
                                return features_v,sheet_obj,False
        
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

        print("Beginning Ring construction.")
        sys.stdout.flush()
        bias_to_TP = (BeNTF2dict["Opening"] == "Tropical")
        bias_to_CH_not_TP = ((not (BeNTF2dict["Opening"] == "Tropical")) and BeNTF2dict["has_cHelix"])

        Ring = CreateRingObjFromDict(ring_data_dict,db=db,bias_TP=bias_to_TP,bias_CH_not_TP=bias_to_CH_not_TP)
        success = AssembleRing(POSE,Ring,ring_trials=25,\
                              prefix=prefix,trial=trial,db=db,options=options)

        if not success: print("Ring construction failed after input attempts") ; return features_v,Ring,False

        can_have_cHelix = NTF2CanHaveCTermH(POSE, ring_obj=Ring)
        tropical_capacity = RingCanBeTP_NTF2(POSE,Ring)

        if BeNTF2dict["Opening"] == "Tropical" and not tropical_capacity:
            print("NTF2 should be able to accomodate alternative opening, but ring is not able, restarting")
            success = False

        if (BeNTF2dict["Opening"] != "Tropical") and BeNTF2dict["has_cHelix"] and (not can_have_cHelix):
            print("NTF2 should be able to accomodate C-helix, but ring is not able, restarting")
            success = False
        
        if not success: print("Ring construction failed because it can't accomodate the C-terminal helix"); return features_v,Ring,False
 
        print("Beginning N-terminal helix construction")
        BeNTF2_obj = CreateBasicNTF2fromDict(BeNTF2dict,NTF2_is_complete=False,db=db)
        success = AssembleBasicBeNTF2(POSE,BeNTF2_obj,trials=25,\
                                      prefix=prefix,trial=trial,db=db,options=options)

        if not success: return features_v,BeNTF2_obj,False

        if BeNTF2_obj.has_cHelix:
                print("Constructing C-helix")
                success = AddCtermHelix(BeNTF2_obj,BasicNTF2_pose=POSE,trials=25,\
                                        prefix=prefix,trial=trial,db=db,options=options)

        if success: return features_v,BeNTF2_obj,True
        else: return features_v,BeNTF2_obj,False

if __name__ == "__main__":
        ########## Option system: ###########
        argparser = argparse.ArgumentParser(description='Create a BeNTF2 backbone.')
        argparser.add_argument('-database', type=str,help='Git repo directory contaninig necessary *.xml and *.wts files')
        argparser.add_argument('-input_pdb', type=str,help='Input PDB file (or just text file) containing a line starting with "BENTF2DICT" followed by a python dictionary json string')
        argparser.add_argument('-nstruct', type=int,help='Max number of structures to output. May output less if some of them fail.')
        argparser.add_argument('-prefix', type=str,help='Prefix to add to output names.')
        args = argparser.parse_args()
        prefix=''
        if args.prefix:
                prefix = args.prefix
        else:
                prefix = os.path.basename(args.input_pdb).split('_')[0]+'_'+randomword(8)
        print("Using prefix: %s"%prefix)

        db = args.database
        general_flags = ['-rama_map %s/Rama_XPG_3level.txt'%db,\
                '-picking_old_max_score 1',\
                '-restore_talaris_behavior',
                '-out:file:pdb_comments',\
                '-mute all']

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
        features_v = {}
        for trial in range(args.nstruct):
                sys.stdout.flush()
                print("Beginning nstruct %d"%trial)
                if os.path.isfile('%s_%d.pdb'%(prefix,trial)):
                        continue
                success =  False
                for attempt in range(20):
                        print("Beginning attempt %d at nstruct %d"%(attempt,trial))
                        sys.stdout.flush()
                        POSE = pose_from_file( filename='%s/input.pdb'%db )
                        if len(features_lines)==1:
                            feat_string = ' '.join(features_lines[0].split()[1:])
                            features_dict = json.loads(feat_string)
                            arch_len_val = features_dict["dist"]
                            features_v,BeNTF2_obj,success = CreateBeNTF2PoseFromDict(BeNTF2dict,POSE,arch_len=arch_len_val)
                        else:
                            features_v,BeNTF2_obj,success = CreateBeNTF2PoseFromDict(BeNTF2dict,POSE)
                        sys.stdout.flush()
                        if success: break
                if not success:
                        print("Unable to make this NTF2 again in the number of iterations provided... sorry...")
                        continue

                # I decided not to add TRP lock and gly positions, those can be automatically detected later:

                POSE.data().clear()

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
