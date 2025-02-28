import os
import argparse


def main(args: argparse.Namespace):
    for pdb_ccd in os.listdir(args.input_folder):
        if pdb_ccd not in ['7M6C_DAO', '7EHR_1PE', '8Z2W_A1L0T', '7EUZ_JD3', '7Y88_GRG', '7Z2S_IHS', '7UY4_SMI',
                           '8ARJ_NRR', '8HUS_7YY', '8Y2S_A1LXW', '7X5N_5M5', '7MCJ_FD7', '7W06_ITN', '8D6M_9SL',
                           '7HM7_EPE', '8IND_6JI', '7ZXS_KBO', '7Q7I_9I8', '7OYJ_3DI', '7Q2B_M6H', '7TXK_LW8',
                           '7SCP_8S6', '8SXY_MJ3', '7P5L_MLT', '7UD8_MW0', '7SS6_BKB', '7U3E_L4O', '9FTC_A1IF6',
                           '8UCC_WFQ', '7AYY_SEQ', '5SCM_GHU', '8QEV_CP', '7VBU_6I4', '7S6V_8II', '8R5B_Y32',
                           '7WKL_CAQ', '8W4S_AJR', '8IVG_DTT', '8C14_T0L', '7F21_FMN', '8PNX_PMP', '7SGH_99W',
                           '8OV7_W3W', '7VDU_65L', '7N7G_FCN', '7NU0_DCL', '8RZH_A1H38', '7ZUY_MTA', '8G4A_YL8',
                           '7Q6U_963', '8BN5_ABU', '7VB8_STL', '7F06_0TL', '8CSC_C5P', '6ZK3_RIB', '7L8N_MLI',
                           '7W2W_BGC', '7THI_PGA', '8YHR_A1LYU', '8EFN_4BW', '8J2X_BLA', '7ZRN_PPV', '8SWS_IM5',
                           '7CQS_GJL', '7VNX_GMP', '7FB7_8NF', '8QML_41K', '7SNB_AMP', '7SCW_GSP', '7R7K_25J',
                           '7OLI_8HG', '7WUY_76N', '8Z4R_A1D7Y', '7VTG_FJF', '7ZO2_DQM', '7V6X_5RJ', '7FEW_CMP',
                           '8CVQ_OYH', '7CQR_GEU', '8TNB_IQM', '8H7N_WYW', '8EWZ_X0I', '7XYK_UMP', '5SHO_JIC',
                           '7ZTL_BCN', '8QMZ_W6O', '7E6A_PG4', '8UOK_X5N', '8EH0_QOX', '8AZN_RVK', '7UM8_NQF',
                           '8JB3_NOS', '7VJT_7IJ', '7UAW_MF6', '8AN6_BTN', '7Z62_IEC', '7S90_ADE', '8KH3_VWF',
                           '7FQ1_FOL', '7V29_5JR', '8RBA_YLW', '7EZQ_0O0', '8RCW_YQF', '8PP8_3IA', '7UZ2_OML',
                           '7PCG_9TY', '7MX6_SNU', '8HFN_XGC', '9GIC_A1ILQ', '7TYP_KUR', '7V36_5P0', '8A1H_DLZ',
                           '7PRI_7TI', '8D5H_HHR', '8ULE_5N6', '8CNC_6D6', '7R9N_F97', '8C0G_SRC', '8CZL_GSM',
                           '8A08_GSH', '8UWN_XQX', '7T3E_SLB', '8B2A_PVI', '8OZ8_NM2', '8PTS_EIF', '8TFO_MEV',
                           '8EGF_WIK', '7NAJ_1LK', '9BFB_A1AN9', '8OI9_PIN', '7P7W_NDG', '7EJ3_DPO', '8R2B_XN0',
                           '8WCH_PCA', '7Y1H_F99', '8BXC_S4R', '8WTQ_KUO', '7EJI_PEG', '7HKD_8BT', '8A46_L5F',
                           '8SQ7_4J6', '8RV7_A1H3E', '7ZF0_DHR', '7W8X_8IK', '7TB0_FFQ', '7WZX_7PK', '7EER_ZCW',
                           '8Q0T_IJJ', '9F93_A1IAV', '8RGK_GOQ', '8GEZ_ZMB', '7ZHN_VP7', '7ULC_56B', '7TBU_S3P',
                           '7OZ8_NGS', '8E9O_MAE', '8FWZ_YCE', '7VYJ_CA0', '8QHG_V8I', '7OF8_VCQ', '8RU8_I74',
                           '7UQ3_O2U', '7ZB3_IVG', '8A7L_ZOL', '8CHL_USV', '8AAU_LH0', '8QFW_UK9', '7POM_7VZ',
                           '7DZ3_FUD', '8PO4_26X', '7YQR_NCA', '8DJU_SIV', '7FD6_4QU', '7QZL_NAP', '7HOC_K3J',
                           '8EPR_NAG', '7ULJ_42C', '7P42_5I8', '7LTY_YD7', '8SZ6_XII', '9BRX_A1ASL', '8HKA_UB7',
                           '8JDB_UIT', '7ZWL_K43', '8BSD_TBN', '8PXH_TAU', '7PRM_81I', '8VKA_A1AB7', '7ROU_66I',
                           '7BBN_GDP', '7MWU_ZPM', '7LCU_XTA', '8I3P_OS0', '7V8W_MLA', '8CB3_MPO', '7QPP_VDX',
                           '7FTG_1SY', '7ELT_TYM', '7TE8_P0T', '7P4C_5OV', '7QGP_DJ8', '5SM7_LPU', '7R1M_BTB',
                           '8CIW_ZKK', '7TSG_KL0']:
            continue
        os.chdir(os.path.join(args.input_folder, pdb_ccd))
        os.system(f"python evaluate.py --pdb_ccd {pdb_ccd}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", type=str, required=True, help="Path to the input folder")
    args = parser.parse_args()

    main(args)
