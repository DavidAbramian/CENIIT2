import argparse
import glob
import os
import random

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from skimage import measure


COLORS = {
    "red" : [255, 0, 0],
    "tomato" : [255, 99, 71],
    "orange" : [255, 165, 0],
    "yellow-green" : [154, 205, 50],
    "lime" : [0, 255, 0],
    "cyan" : [0, 255, 255],
    "turquoise" : [64, 224, 208],
    "deep-blue-sky" : [0, 191, 255],
    "blue" : [0, 0, 255],
    "dark-violet" : [148, 0, 211],
    "magenta" : [255, 0, 255],
    "brown" : [165, 42, 42],
    "tan" : [210, 160, 120],
    "pink" : [255, 192, 203],
}


def convert(subject: str, str_stats: str, str_run: str, series_number: int, contour_level: float):

    # subject = "17904"
    # str1 = "z5_c0"
    # str2 = "4mm_m"

    # subject = "18582"
    # str_stats = "z5_c0"
    # str_run = "4mm_m"

    activity_names = [f"{x}_{str_run}" for x in  ["tumor", "fing", "foot", "lips", "verb"]]
    colors = ["red", "lime", "orange", "cyan", "magenta"]
    ROIType = ["PTV", "AVOIDANCE", "AVOIDANCE", "AVOIDANCE", "AVOIDANCE"]
    contourLevelStr = f"{contour_level*10:02.0f}"

    # Tumor mask
    tumor_nifti_file = os.path.join(subject, "anat", "tumor_mask_new.nii.gz")

    # Activation maps
    activations_nifti_path = os.path.join(subject, "fmri_stats", str_stats)
    activations_nifti_files = sorted(glob.glob(os.path.join(activations_nifti_path, f"*{str_run}.nii.gz")))

    nifti_files = [tumor_nifti_file] + activations_nifti_files
    n_activations_files = len(nifti_files)

    # T1w DICOM files
    t1_dicom_path = os.path.join(subject, "dicom", f"{subject}_dicom")
    t1_dicom_files = sorted(glob.glob(os.path.join(t1_dicom_path, "*.dcm")))
    n_dicom_files = len(t1_dicom_files)

    activity_names = activity_names[:n_activations_files]
    colors = colors[:n_activations_files]
    ROIType = ROIType[:n_activations_files]

    # Output
    output_rstruct_path = os.path.join(subject, "dicom", f"{subject}_dicom_RTSTRUCT")
    if not os.path.exists(output_rstruct_path):
        os.mkdir(output_rstruct_path)
    output_file_name = f"{subject}_{str_run}_{contourLevelStr}"


    #---------------
    # First DICOM part
    #---------------

    dicom_params_dict = get_dicom_params(t1_dicom_files)
    dicom_params_dict["activity_names"] = activity_names
    dicom_params_dict["colors"] = colors
    dicom_params_dict["ROIType"] = ROIType
    dicom_params_dict["series_number"] = series_number
    dicom_params_dict["series_number"] = series_number
    dicom_params_dict["output_file_name"] = output_file_name
    dicom_params_dict["str_stats"] = str_stats
    dicom_params_dict["str_run"] = str_run

    #---------------
    # NIFTI part
    #---------------

    # Find position of first slice
    patientPosition = np.array(dicom_params_dict["ImagePositionPatient"][0], ndmin=2)
    patientPosition[0,2] = np.min([a[2] for a in dicom_params_dict["ImagePositionPatient"]])

    # List to store contours
    all_file_ROI_list = []
    n_ROIs_per_file = np.zeros(n_activations_files, dtype=int)

    # Find ROI contours
    for f, nifti_file in enumerate(nifti_files):

        file_ROIs = extract_contours(nifti_file, patientPosition, contour_level)
        all_file_ROI_list.append(file_ROIs)
        n_ROIs_per_file[f] = len(file_ROIs)

        print(f"File {nifti_file} has {n_ROIs_per_file[f]} ROIs")

    # Plot ROI contours
    if False:

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for ROI in file_ROIs:

            col = np.random.rand(3)
            for _, contour in ROI:

                # ax.plot(contour[0::3], contour[1::3], contour[2::3],'-',color=np.array(list(colors.values()))[r] / 255)
                ax.plot(contour[0::3], contour[1::3], contour[2::3],'-',color=col)

    # Plot ROI contours
    if False:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for f in range(n_activations_files):

            for ROI in all_file_ROI_list[f]:

                for _, contour in ROI:

                    ax.plot(contour[0::3], contour[1::3], contour[2::3],'-',color=colors[f])

    #---------------
    # Second DICOM part (RTstruct)
    #---------------

    ds_rt = create_rtstruct_file(dicom_params_dict, all_file_ROI_list)

    # Save RTSTRUCT file
    ds_rt.save_as(os.path.join(output_rstruct_path, f"{output_file_name}_RTS.dcm"))


def get_dicom_params(t1_dicom_files):

    dicom_dict = dict()
    dicom_dict["FileName"] = []
    dicom_dict["SOPInstanceUID"] = []
    dicom_dict["ImagePositionPatient"] = []

    for file in range(len(t1_dicom_files)):

        ds_t1 = pydicom.dcmread(t1_dicom_files[file], stop_before_pixels=True)

        dicom_dict["FileName"].append(t1_dicom_files[file])
        dicom_dict["SOPInstanceUID"].append(ds_t1.SOPInstanceUID)
        dicom_dict["ImagePositionPatient"].append(ds_t1.ImagePositionPatient)

    dicom_dict["NFiles"] = len(t1_dicom_files)
    dicom_dict["SOPClassUID"] = ds_t1.SOPClassUID
    dicom_dict["StudyInstanceUID"] = ds_t1.StudyInstanceUID
    dicom_dict["StudyDate"] = ds_t1.StudyDate
    dicom_dict["StudyTime"] = ds_t1.StudyTime
    dicom_dict["PatientName"] = ds_t1.PatientName
    dicom_dict["PatientID"] = ds_t1.PatientID
    dicom_dict["PatientBirthDate"] = ds_t1.PatientBirthDate
    dicom_dict["PatientSex"] = ds_t1.PatientSex
    dicom_dict["FrameOfReferenceUID"] = ds_t1.FrameOfReferenceUID
    dicom_dict["SeriesInstanceUID"] = ds_t1.SeriesInstanceUID

    return dicom_dict


def extract_contours(nifti_file, patient_position, contour_level):

    # Load nifti volume
    nii = nib.load(nifti_file)
    volume = nii.get_fdata().astype(float)

    if not len(volume.shape)==3:
        raise Exception('Dimension not supported.')

    pixel_size = np.array(nii.header.get_zooms(), ndmin=2)

    # Label connected components
    volume_labels = measure.label(volume, connectivity=1)
    n_ROIs = volume_labels.max()

    file_ROIs = []
    for ROI in range(1, n_ROIs+1):

        volume_ROI = volume_labels == ROI

        ROI_contours = []
        # Loop over slices in volume, get contours for each slice
        for z in range(volume.shape[2]):

            image = np.fliplr(volume_ROI[:,:,z])

            if not np.any(image):
                continue

            # Get contours in this slice using scikit-image
            contours = measure.find_contours(image, contour_level)

            # Save contours for later use
            for contour in contours:

                # contour = reduce_contour(contour)
                # contour = augment_contour(contour, 5)

                # Remove last (repeated) point from contour
                contour = contour[:-1,:]

                # Add z coordinate
                contour = np.append(contour, z * np.ones((contour.shape[0],1)), axis=1)

                # Transform coordinates from voxels to mm
                contour = contour * pixel_size + patient_position

                # Store unrolled coordinates
                coordinates = contour.reshape(-1)
                ROI_contours.append([z, coordinates])

        file_ROIs.append(ROI_contours)

    return file_ROIs


def reduce_contour(contour):

    # c = np.array([[0, 0], [0,1], [0,2], [0,3], [1,3], [2,3], [3,3], [3,2], [3,1], [3,0], [0,0]])
    c = contour
    d = np.diff(c, axis=0)
    dd = np.abs(np.diff(d, axis=0)).sum(axis=1)
    c2 = np.concatenate((c[[0],:], c[1:-1][dd != 0, :], c[[-1],:]), axis=0)

    # plt.figure(1)
    # plt.clf()
    # plt.subplot(1,2,1)
    # plt.plot(c[:,0], c[:,1], '.-')

    # plt.subplot(1,2,2)
    # plt.plot(c2[:,0], c2[:,1], '.-')
    # plt.show()

    # plt.figure(1)
    # plt.clf()
    # plt.plot(c[:,0], c[:,1], '.-')
    # plt.plot(c2[:,0], c2[:,1], '.-')
    # plt.show()

    return c2


def augment_contour(contour, factor):

    # c = np.array([[0, 0], [0,1], [0,2], [0,3], [1,3], [2,3], [3,3], [3,2], [3,1], [3,0], [0,0]])
    c = contour
    d = np.diff(c, axis=0)

    c2 = np.zeros((factor*c.shape[0]-(factor-1), 2))
    c2[::factor] = c

    for i in range(1, factor+1):
        c2[i::factor] = c[:-1] + d * i/factor

    # plt.figure(1)
    # plt.clf()
    # plt.subplot(1,2,1)
    # plt.plot(c[:,0], c[:,1], '.-')

    # plt.subplot(1,2,2)
    # plt.plot(c2[:,0], c2[:,1], '.-')
    # plt.show()

    # plt.figure(1)
    # plt.clf()
    # plt.plot(c[:,0], c[:,1], '.-')
    # plt.plot(c2[:,0], c2[:,1], '.-')
    # plt.show()

    return c2


def create_rtstruct_file(dicom_params_dict, all_file_ROI_list):

    n_dicom_files = dicom_params_dict["NFiles"]
    n_activations_files = len(dicom_params_dict["ROIType"])
    n_ROIs_per_file = [len(a) for a in all_file_ROI_list]

    ds_rt = pydicom.Dataset()
    ds_rt.is_implicit_VR = False
    ds_rt.is_little_endian = True

    # SOP Common Module Attributes
    ds_rt.SpecificCharacterSet = "ISO_IR 192"  # (0008,0005), Type 1C
    ds_rt.InstanceCreationDate = dicom_params_dict["StudyDate"]  # (0008, 0012) Type 3
    ds_rt.InstanceCreationTime = dicom_params_dict["StudyTime"]  # (0008, 0013) Type 3
    ds_rt.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'  # (0008, 0016) Type 1, RTSTRUCT
    ds_rt.SOPInstanceUID = f"1.2.826.0.1.123123123.{random.randint(0,99999999)}.1.{random.randint(0,99999999)}"  # (0008,0018) Type 1

    # General Study Module Attributes
    ds_rt.StudyInstanceUID = dicom_params_dict["StudyInstanceUID"] # (0020,000D), Type 1
    ds_rt.StudyDate = dicom_params_dict["StudyDate"]  # (0008,0020), Type 2
    ds_rt.StudyTime = dicom_params_dict["StudyTime"]  # (0008,0030), Type 2
    ds_rt.ReferringPhysicianName = "refPhys"  # (0008, 0090), Type 2
    ds_rt.StudyID = "studyId"  # (0020,0010), Type 2
    ds_rt.AccessionNumber = ""  # (0008,0050), Type 2
    ds_rt.StudyDescription = "studyDesc"  # (0008,1030), Type 3

    # General Series Module Attributes
    ds_rt.Modality = 'RTSTRUCT'  # (0008,0060), Type 1
    ds_rt.SeriesInstanceUID = f"1.2.826.0.1.3680043.2.1125.{random.randint(0,99999999)}.1.{random.randint(0,99999999)}"  # (0020,000E) Type 1
    ds_rt.SeriesNumber = dicom_params_dict["series_number"] # (0020,0011), Type 2
    ds_rt.SeriesDescription = dicom_params_dict["output_file_name"]  # (0008.103E), Type 3

    # General Equipement Module Attributes
    ds_rt.Manufacturer = "manuf"  # (0008,0070), Type 2
    # ds_rt.InstitutionName = ""  # (0008,0080), Type 3
    # ds_rt.InstitutionAddress = ""  # (0008,0081), Type 3
    # ds_rt.StationName = ""  # (0008,1010), Type 3
    # ds_rt.InstitutionalDepartmentName = ""  # (0008,1040), Type 3
    # ds_rt.ManufacturerModelName = ""  # (0008,1090), Type 3
    # ds_rt.DeviceSerialNumber = ""  # (0018,1000), Type 3
    # ds_rt.SoftwareVersion = ""  # (0018,1020), Type 3

    # Patient Module Attributes
    ds_rt.PatientName = dicom_params_dict["PatientName"]  # (0010,0010), Type 2
    ds_rt.PatientID = dicom_params_dict["PatientID"]  # (0010,0020), Type 2
    ds_rt.PatientBirthDate = dicom_params_dict["PatientBirthDate"]  # (0010,0030), Type 2
    ds_rt.PatientSex = dicom_params_dict["PatientBirthDate"]  # (0010,0040), Type 2

    # Frame of Reference Module Attributes
    ds_rt.FrameOfReferenceUID = dicom_params_dict["FrameOfReferenceUID"]  # (0020, 0052), Type 1
    ds_rt.PositionReferenceIndicator = ""  # (0020, 0052), Type 2

    # Structure Set Module Attributes
    ds_rt.StructureSetLabel = dicom_params_dict["str_run"]  # (3006,0002), Type 1
    ds_rt.StructureSetDate = ""  # (3006,0008), Type 2
    ds_rt.StructureSetTime = ""  # (3006,0009), Type 2


    # Referenced Frame of Reference Sequence
    refd_frame_of_ref_sequence = Sequence()
    ds_rt.ReferencedFrameOfReferenceSequence = refd_frame_of_ref_sequence

    # Referenced Frame of Reference Sequence: Referenced Frame of Reference 1
    refd_frame_of_ref1 = Dataset()
    refd_frame_of_ref1.FrameOfReferenceUID = dicom_params_dict["FrameOfReferenceUID"]

    # RT Referenced Study Sequence
    rt_refd_study_sequence = Sequence()
    refd_frame_of_ref1.RTReferencedStudySequence = rt_refd_study_sequence

    # RT Referenced Study Sequence: RT Referenced Study 1
    rt_refd_study1 = Dataset()
    rt_refd_study1.ReferencedSOPClassUID = dicom_params_dict["SOPClassUID"]
    rt_refd_study1.ReferencedSOPInstanceUID = dicom_params_dict["StudyInstanceUID"]

    # RT Referenced Series Sequence
    rt_refd_series_sequence = Sequence()
    rt_refd_study1.RTReferencedSeriesSequence = rt_refd_series_sequence

    # RT Referenced Series Sequence: RT Referenced Series 1
    rt_refd_series1 = Dataset()
    rt_refd_series1.SeriesInstanceUID = dicom_params_dict["SeriesInstanceUID"]

    # Contour Image Sequence
    contour_image_sequence = Sequence()
    rt_refd_series1.ContourImageSequence = contour_image_sequence

    for file in range(n_dicom_files):
        # Contour Image Sequence: Contour Image
        contour_image = Dataset()
        contour_image.ReferencedSOPClassUID = dicom_params_dict["SOPClassUID"]
        contour_image.ReferencedSOPInstanceUID = dicom_params_dict["SOPInstanceUID"][file]
        contour_image_sequence.append(contour_image)

    rt_refd_series_sequence.append(rt_refd_series1)
    rt_refd_study_sequence.append(rt_refd_study1)
    refd_frame_of_ref_sequence.append(refd_frame_of_ref1)

    # Structure Set ROI Sequence
    structure_set_roi_sequence = Sequence()
    ds_rt.StructureSetROISequence = structure_set_roi_sequence

    ROI_counter = 0
    for file in range(n_activations_files):
        for ROI in range(1, n_ROIs_per_file[file]+1):
            ROI_counter += 1

            # Structure Set ROI Sequence: Structure Set ROI
            structure_set_roi = Dataset()
            structure_set_roi.ROINumber = str(ROI_counter)
            structure_set_roi.ReferencedFrameOfReferenceUID = dicom_params_dict["FrameOfReferenceUID"]
            structure_set_roi.ROIName = f"{dicom_params_dict['activity_names'][file]}_{ROI:03}"
            structure_set_roi.ROIGenerationAlgorithm = ""
            structure_set_roi_sequence.append(structure_set_roi)


    # ROI Contour module
    # ROI Contour Sequence
    roi_contour_sequence = Sequence()
    ds_rt.ROIContourSequence = roi_contour_sequence

    ROI_counter = 0
    for file in range(n_activations_files):

        # Loop over ROI contour sequences
        for ROI in range(n_ROIs_per_file[file]):
            ROI_counter += 1

            # ROI Contour Sequence: ROI Contour 1
            roi_contour = Dataset()
            roi_contour.ReferencedROINumber = ROI_counter
            roi_contour.ROIDisplayColor = COLORS[dicom_params_dict["colors"][file]]
            roi_contour_sequence.append(roi_contour)

            # Contour Sequence
            contour_sequence = Sequence()
            roi_contour.ContourSequence = contour_sequence

            for z, contour_coordinates in all_file_ROI_list[file][ROI]:

                # Contour Sequence: Contour 1
                contour = Dataset()

                # Contour Image Sequence
                contour_image_sequence = Sequence()
                contour.ContourImageSequence = contour_image_sequence

                # Contour Image Sequence: Contour Image 1
                contour_image = Dataset()
                contour_image.ReferencedSOPClassUID = dicom_params_dict["SOPClassUID"]  # '1.2.840.10008.5.1.4.1.1.2'
                contour_image.ReferencedSOPInstanceUID = dicom_params_dict["SOPInstanceUID"][z] # '1.3.6.1.4.1.9590.100.1.2.76071554513024464020636223132290799275'
                contour_image_sequence.append(contour_image)

                contour.ContourGeometricType = 'CLOSED_PLANAR'
                contour.NumberOfContourPoints = len(contour_coordinates) / 3
                contour.ContourData = contour_coordinates.tolist()
                contour_sequence.append(contour)


    # RT ROI Observations Module
    # RT ROI Observations Sequence
    rtroi_observations_sequence = Sequence()
    ds_rt.RTROIObservationsSequence = rtroi_observations_sequence

    ROI_counter = 0
    for file in range(n_activations_files):
        for ROI in range(n_ROIs_per_file[file]):
            ROI_counter += 1

            # RT ROI Observations Sequence: RT ROI Observations 1
            rtroi_observations = Dataset()
            rtroi_observations.ObservationNumber = str(ROI_counter)
            rtroi_observations.ReferencedROINumber = str(ROI_counter)
            rtroi_observations.ROIObservationLabel = ''
            rtroi_observations.RTROIInterpretedType = dicom_params_dict["ROIType"][file]
            rtroi_observations.ROIInterpreter = ''
            rtroi_observations_sequence.append(rtroi_observations)

    return ds_rt


def get_parser():
    parser = argparse.ArgumentParser(description='Convert nifti segmentation volume to DICOM RTstruct')
    # Positional arguments.
    # parser.add_argument("input_nifti", help="Path to input NIFTI image containing segmentation mask")
    # parser.add_argument("input_dicom", help="Path to original DICOM files (used as template)")
    # parser.add_argument("output_dicom", help="Path to output DICOM RTSTRUCT image")
    # parser.add_argument("ROI_name", help="Name of this ROI")
    # parser.add_argument("ROI_color", help="Color of this ROI")

    parser.add_argument("subject", help="Subject ID")
    parser.add_argument("str1", help="fMRI run string, e.g. z5_c0")
    parser.add_argument("str2", help="activation map string, e.g. 6mm_m4")
    parser.add_argument("series_number", type=int, help="series number")
    parser.add_argument("-l", "--level", metavar="", type=float, default=0.5, help="threshold for finding contours [0-1], default 0.5")

    return parser.parse_args()

if __name__ == "__main__":
    p = get_parser()
    convert(p.subject, p.str1, p.str2, p.series_number, p.level)
