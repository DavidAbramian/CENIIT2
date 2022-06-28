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

def find_first_slice_position(dcms):

    patientStartingZ = 0
    for idx, dcm in enumerate(dcms):
        ds = pydicom.dcmread(dcm, stop_before_pixels=True)
        if not 'ImagePositionPatient' in ds or ds.ImagePositionPatient is None:
            continue
        if ds.ImagePositionPatient[2] <= patientStartingZ or idx==0:
            patientStartingZ = ds.ImagePositionPatient[2]

    return patientStartingZ


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


# def convert(input_nifti_path: str, input_dicom_path: str, output_dicom_path: str, ROI_name: str, ROI_color: str):
def convert(subject: str, str1: str, str2: str):

    # subject = "17904"
    # str1 = "z5_c0"
    # str2 = "4mm_m"

    subject = "18582"
    str1 = "z5_c0"
    str2 = "4mm_m"

    input_dicom_path = os.path.join(subject, "dicom", f"{subject}_dicom")
    output_dicom_path = os.path.join(subject, "dicom", f"{subject}_dicom_RTSTRUCT")
    # out_name = f"{subject}_{str2}"
    out_name = f"{subject}_{str2}_04_6con_2"

    nifti_path = os.path.join(subject, "fmri_stats", str1)
    nifti_files = sorted(glob.glob(os.path.join(nifti_path, f"*{str2}.nii.gz")))
    n_files = len(nifti_files)

    #---------------
    # First DICOM part
    #---------------

    # Get number of DICOM files in DICOM path
    dicom_files = next(os.walk(input_dicom_path))[2]
    n_dicom_files = len(dicom_files)
    # n_ROIs = 1   # The whole volume is 1 ROI, assuming 1 tumour per patient

    # Load template DICOM file header (first file)
    ds = pydicom.dcmread(os.path.join(input_dicom_path, dicom_files[0]), stop_before_pixels=True)

    x_pixel_size = ds.PixelSpacing[0]
    y_pixel_size = ds.PixelSpacing[1]
    z_pixel_size = ds.SliceThickness

    # print(f"Voxel size: {x_pixel_size} x {y_pixel_size} x {z_pixel_size} mm")

    # Find position of first slice
    patientPosition = ds.ImagePositionPatient
    patientStartingZ = find_first_slice_position([os.path.join(input_dicom_path, x) for x in dicom_files])

    # print(f"Patient position: {patientPosition[:2]}")
    # print(f"First slice at {patientStartingZ}")

    #---------------
    # NIFTI part
    #---------------

    all_file_ROI_list = []
    n_ROIs_per_file = np.zeros(n_files, dtype=int)
    activity_names = [f"{x}_{str2}" for x in  ["fing", "foot", "lips", "verb"]]
    colors = ["lime", "orange", "cyan", "magenta"]

    for f, nifti_file in enumerate(nifti_files):

        # Load nifti volume
        nii = nib.load(nifti_file)
        volume = nii.get_fdata().astype(float)

        if not len(volume.shape)==3:
            raise Exception('Dimension not supported.')

        # Label connected components
        volume_labels = measure.label(volume, connectivity=1)
        n_ROIs = volume_labels.max()
        n_ROIs_per_file[f] = n_ROIs

        file_ROIs = []

        for ROI in range(1, n_ROIs+1):

            volume_ROI = volume_labels == ROI

            ROI_slices = []

            # Loop over slices in volume, get contours for each slice
            for z in range(volume.shape[2]):

                image = volume_ROI[:,:,z]
                image = np.fliplr(image)

                if not np.any(image):
                    ROI_slices.append([])
                else:
                    slice_contours = []

                    # Get contours in this slice using scikit-image
                    # contours = measure.find_contours(image, 0.5)
                    # contours = measure.find_contours(image, 0.6)
                    contours = measure.find_contours(image, 0.4)

                    # Save contours for later use
                    for contour in contours:

                        # contour = reduce_contour(contour)
                        # contour = augment_contour(contour, 5)


                        # Remove last (repeated) point from contour
                        # contour = contour[:-1,:]

                        # Add patient position offset
                        n_coordinates = contour.shape[0]
                        z_coordinates = z * np.ones((n_coordinates,1))
                        contour = np.append(contour, z_coordinates, -1)

                        # Storing coordinates as mm instead of as voxels
                        contour[:,0] = contour[:,0] * x_pixel_size + patientPosition[0]
                        contour[:,1] = contour[:,1] * y_pixel_size + patientPosition[1]
                        contour[:,2] = contour[:,2] * z_pixel_size + patientStartingZ

                        # Unroll coordinates
                        coordinates = contour.reshape(-1)
                        slice_contours.append(coordinates)

                    ROI_slices.append(slice_contours)

            file_ROIs.append(ROI_slices)

        all_file_ROI_list.append(file_ROIs)

        # # Print number of slices in each ROI
        # nSlicesPerROI = [np.flatnonzero([len(y) for y in x]).shape[0] for x in file_ROIs]
        # for i, x in enumerate(nSlicesPerROI):
        #     print(f"{i+1}: {x}")

        print(f"File {nifti_file} has {len(file_ROIs)} ROIs")

    #print("All coordinates has length ",len(AllCoordinates))
    #print("All coordinates slice 0 has length ",len(AllCoordinates[0]))
    #print("All coordinates slice 1 has length ",len(AllCoordinates[1]))
    #print("All coordinates slice 1 contour 1 has length ",len(AllCoordinates[1][1]))
    #print("Coordinates are ",AllCoordinates[1][1])

    # Plot contours
    if False:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for ROI, ROI in enumerate(file_ROIs):

            # if r != 11:
            #     continue

            col = color=np.random.rand(3)
            for slice in ROI:
                for contour in slice:
                    # ax.plot(contour[0::3], contour[1::3], contour[2::3],'-',color=np.array(list(colors.values()))[r] / 255)
                    ax.plot(contour[0::3], contour[1::3], contour[2::3],'-',color=col)


    #---------------
    # Second DICOM part (RTstruct)
    #---------------

    # Referenced Frame of Reference Sequence
    refd_frame_of_ref_sequence = Sequence()
    ds.ReferencedFrameOfReferenceSequence = refd_frame_of_ref_sequence

    # Referenced Frame of Reference Sequence: Referenced Frame of Reference 1
    refd_frame_of_ref1 = Dataset()
    refd_frame_of_ref1.FrameOfReferenceUID = ds.FrameOfReferenceUID # '1.3.6.1.4.1.9590.100.1.2.138467792711241923028335441031194506417'

    # RT Referenced Study Sequence
    rt_refd_study_sequence = Sequence()
    refd_frame_of_ref1.RTReferencedStudySequence = rt_refd_study_sequence

    # RT Referenced Study Sequence: RT Referenced Study 1
    rt_refd_study1 = Dataset()
    rt_refd_study1.ReferencedSOPClassUID = ds.SOPClassUID # '1.2.840.10008.5.1.4.1.1.481.3'
    rt_refd_study1.ReferencedSOPInstanceUID = ds.SOPInstanceUID # '1.3.6.1.4.1.9590.100.1.2.201285932711485367426568006803977990318'

    # RT Referenced Series Sequence
    rt_refd_series_sequence = Sequence()
    rt_refd_study1.RTReferencedSeriesSequence = rt_refd_series_sequence

    # RT Referenced Series Sequence: RT Referenced Series 1
    rt_refd_series1 = Dataset()
    rt_refd_series1.SeriesInstanceUID = ds.SeriesInstanceUID  # '1.3.6.1.4.1.9590.100.1.2.170217758912108379426621313680109428629'

    # Contour Image Sequence
    contour_image_sequence = Sequence()
    rt_refd_series1.ContourImageSequence = contour_image_sequence

    # Loop over all DICOM images
    for image in range(1, n_dicom_files+1):
        dstemp = pydicom.dcmread(os.path.join(input_dicom_path, dicom_files[image-1]),stop_before_pixels=True)
        # Contour Image Sequence: Contour Image
        contour_image = Dataset()
        contour_image.ReferencedSOPClassUID = dstemp.SOPClassUID      # '1.2.840.10008.5.1.4.1.1.2'
        contour_image.ReferencedSOPInstanceUID = dstemp.SOPInstanceUID # '1.3.6.1.4.1.9590.100.1.2.257233736012685791123157667031991108836'
        contour_image_sequence.append(contour_image)

    rt_refd_series_sequence.append(rt_refd_series1)
    rt_refd_study_sequence.append(rt_refd_study1)
    refd_frame_of_ref_sequence.append(refd_frame_of_ref1)

    # Structure Set ROI Sequence
    structure_set_roi_sequence = Sequence()
    ds.StructureSetROISequence = structure_set_roi_sequence

    # Loop over ROIs
    ROI_counter = 0
    for file in range(n_files):
        for ROI in range(1, n_ROIs_per_file[file]+1):
            ROI_counter += 1

            # Structure Set ROI Sequence: Structure Set ROI
            structure_set_roi = Dataset()
            structure_set_roi.ROINumber = str(ROI_counter)
            structure_set_roi.ReferencedFrameOfReferenceUID = ds.FrameOfReferenceUID # '1.3.6.1.4.1.9590.100.1.2.138467792711241923028335441031194506417'
            structure_set_roi.ROIName = f"{activity_names[file]}_{ROI:03}"
            structure_set_roi.ROIGenerationAlgorithm = 'AndersNicePythonScript'
            structure_set_roi_sequence.append(structure_set_roi)

    # ROI Contour Sequence
    roi_contour_sequence = Sequence()
    ds.ROIContourSequence = roi_contour_sequence

    ROI_counter = 0
    for file in range(n_files):

        # Loop over ROI contour sequences
        for ROI in range(1,n_ROIs_per_file[file]+1):
            ROI_counter += 1

            # ROI Contour Sequence: ROI Contour 1
            roi_contour = Dataset()
            roi_contour.ROIDisplayColor = COLORS[colors[file]]

            # Contour Sequence
            contour_sequence = Sequence()
            roi_contour.ContourSequence = contour_sequence

            # Loop over slices in volume (ROI)
            for z in range(volume.shape[2]):

                # Should Contour Sequence be inside this loop?
                #contour_sequence = Sequence()
                #roi_contour.ContourSequence = contour_sequence

                # Loop over contour sequences in this slice
                numberOfContoursInThisSlice = len(all_file_ROI_list[file][ROI-1][z])
                for c in range(numberOfContoursInThisSlice):

                    currentCoordinates = all_file_ROI_list[file][ROI-1][z][c]

                    # Contour Sequence: Contour 1
                    contour = Dataset()

                    # Contour Image Sequence
                    contour_image_sequence = Sequence()
                    contour.ContourImageSequence = contour_image_sequence

                    # Load the corresponding dicom file to get the SOPInstanceUID
                    dstemp = pydicom.dcmread(os.path.join(input_dicom_path, dicom_files[z]), stop_before_pixels=True)

                    # Contour Image Sequence: Contour Image 1
                    contour_image = Dataset()
                    contour_image.ReferencedSOPClassUID = dstemp.SOPClassUID  # '1.2.840.10008.5.1.4.1.1.2'
                    contour_image.ReferencedSOPInstanceUID = dstemp.SOPInstanceUID # '1.3.6.1.4.1.9590.100.1.2.76071554513024464020636223132290799275'
                    contour_image_sequence.append(contour_image)

                    contour.ContourGeometricType = 'CLOSED_PLANAR'
                    contour.NumberOfContourPoints = len(currentCoordinates)
                    contour.ContourData = currentCoordinates.tolist()
                    contour_sequence.append(contour)

            roi_contour.ReferencedROINumber = ROI_counter
            roi_contour_sequence.append(roi_contour)


    # RT ROI Observations Sequence
    rtroi_observations_sequence = Sequence()
    ds.RTROIObservationsSequence = rtroi_observations_sequence

    # Loop over ROI observations
    ROI_counter = 0
    for file in range(n_files):
        for ROI in range(1, n_ROIs_per_file[file]+1):
            ROI_counter += 1

            # RT ROI Observations Sequence: RT ROI Observations 1
            rtroi_observations = Dataset()
            rtroi_observations.ObservationNumber = str(ROI_counter)
            rtroi_observations.ReferencedROINumber = str(ROI_counter)
            rtroi_observations.ROIObservationLabel = ''
            rtroi_observations.RTROIInterpretedType = 'AVOIDANCE'
            rtroi_observations.ROIInterpreter = ''
            rtroi_observations_sequence.append(rtroi_observations)

    # Add RTSTRUCT specifics
    ds.Modality = 'RTSTRUCT' # So the GammaPlan software can recognize RTSTRUCT
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3' # So the GammaPlan software can recognize RTSTRUCT

    random_str_1 = "%0.8d" % random.randint(0,99999999)
    random_str_2 = "%0.8d" % random.randint(0,99999999)
    ds.SeriesInstanceUID = "1.2.826.0.1.3680043.2.1125."+random_str_1+".1"+random_str_2 # Just some random UID

    # Rename series for identification in GammaPlan
    ds.SeriesDescription = out_name

    # Save RTSTRUCT file
    ds.save_as(os.path.join(output_dicom_path, f"{out_name}_RTSTRUCT.dcm"))

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

    return parser.parse_args()

if __name__ == "__main__":
    p = get_parser()
    convert(p.subject, p.str1, p.str2)
