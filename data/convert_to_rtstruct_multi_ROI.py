import argparse
import os
import random

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from skimage import measure

def concatenate_coordinates(coordinates_x, coordinates_y, coordinates_z):

    vector = np.zeros(coordinates_x.shape[0]*3)
    vector[0::3] = coordinates_x
    vector[1::3] = coordinates_y
    vector[2::3] = coordinates_z

    return vector

def find_first_slice_position(dcms):

    patientStartingZ = 0
    for idx, dcm in enumerate(dcms):
        ds = pydicom.dcmread(dcm, stop_before_pixels=True)
        if not 'ImagePositionPatient' in ds or ds.ImagePositionPatient is None:
            continue
        if ds.ImagePositionPatient[2] <= patientStartingZ or idx==0:
            patientStartingZ = ds.ImagePositionPatient[2]

    return patientStartingZ

def convert(input_nifti_path: str, input_dicom_path: str, output_dicom_path: str, ROI_name: str, ROI_color: str):

    colors = {
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

    #---------------
    # First DICOM part
    #---------------

    # Get number of DICOM files in DICOM path
    dicomFiles = next(os.walk(input_dicom_path))[2]
    numberOfDicomImages = len(dicomFiles)
    numberOfROIs = 1   # The whole volume is 1 ROI, assuming 1 tumour per patient

    # Load template DICOM file header (first file)
    ds = pydicom.dcmread(input_dicom_path + "%s"%dicomFiles[0],stop_before_pixels=True)

    xPixelSize = ds.PixelSpacing[0]
    yPixelSize = ds.PixelSpacing[1]
    zPixelSize = ds.SliceThickness

    print("Each voxel is ",xPixelSize," x ",yPixelSize," x ",zPixelSize)

    # Find position of first slice
    patientPosition = ds.ImagePositionPatient
    patientStartingZ = find_first_slice_position([input_dicom_path + '%s'%_ for _ in dicomFiles])

    print('Patient position is ', patientPosition[:2])
    print('First slice at ', patientStartingZ)

    #---------------
    # NIFTI part
    #---------------

    # Load nifti volume
    nii = nib.load(input_nifti_path)
    volume = nii.get_fdata().astype(float)

    if len(volume.shape)==4:
        volume = volume[...,0]
        print('Assuming the first channel of the input nifti is the seg mask.')
    elif len(volume.shape)==3:
        print('Segmentation mask is same size of the patient image volume.')
    else:
        print('Dimension not supported.')

    # Label connected components
    volumeLabels = measure.label(volume)
    numberOfROIs = volumeLabels.max()

    allContourCoordinates = []

    for r in range(1, numberOfROIs+1):

        volumeROI = volumeLabels == r

        allROIContourCoordinates = []

        # Loop over slices in volume, get contours for each slice
        for z in range(volume.shape[2]):

            contourCoordinatesThisSlice = []

            image = volumeROI[:,:,z]
            image = np.fliplr(image)

            # Get contours in this slice using scikit-image
            # contours = measure.find_contours(image, 0.5)
            contours = measure.find_contours(image, 0.4)

            # Save contours for later use
            for contour in contours:

                contour = contour[:-1,:]

                nCoordinates = contour.shape[0]
                zcoordinates = z * np.ones((nCoordinates,1))

                # Add patient position offset
                reg_contour = np.append(contour, zcoordinates, -1)

                # Assume no other orientations for simplicity
                reg_contour[:,0] = reg_contour[:,0] * xPixelSize + patientPosition[0]
                reg_contour[:,1] = reg_contour[:,1] * yPixelSize + patientPosition[1]
                reg_contour[:,2] = reg_contour[:,2] * zPixelSize + patientStartingZ

                # Storing coordinates as mm instead of as voxels
                coordinates = concatenate_coordinates(*reg_contour.T)

                contourCoordinatesThisSlice.append(coordinates)

            allROIContourCoordinates.append(contourCoordinatesThisSlice)

        allContourCoordinates.append(allROIContourCoordinates)

    #print("All coordinates has length ",len(AllCoordinates))
    #print("All coordinates slice 0 has length ",len(AllCoordinates[0]))
    #print("All coordinates slice 1 has length ",len(AllCoordinates[1]))
    #print("All coordinates slice 1 contour 1 has length ",len(AllCoordinates[1][1]))
    #print("Coordinates are ",AllCoordinates[1][1])

    # Plot contours
    if True:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for r, ROI in enumerate(allContourCoordinates):

            # if r != 11:
            #     continue

            col = color=np.random.rand(3)
            for slice in ROI:
                for contour in slice:
                    # ax.plot(contour[0::3], contour[1::3], contour[2::3],'-',color=np.array(list(colors.values()))[2] / 255)
                    ax.plot(contour[0::3], contour[1::3], contour[2::3],'.-',color=col)


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
    for image in range(1,numberOfDicomImages+1):
        dstemp = pydicom.dcmread(input_dicom_path + "%s"%dicomFiles[image-1],stop_before_pixels=True)
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
    for ROI in range(1,numberOfROIs+1):
        # Structure Set ROI Sequence: Structure Set ROI
        structure_set_roi = Dataset()
        structure_set_roi.ROINumber = str(ROI)
        structure_set_roi.ReferencedFrameOfReferenceUID = ds.FrameOfReferenceUID # '1.3.6.1.4.1.9590.100.1.2.138467792711241923028335441031194506417'
        structure_set_roi.ROIName = f"{ROI_name}_{ROI}"
        structure_set_roi.ROIGenerationAlgorithm = 'AndersNicePythonScript'
        structure_set_roi_sequence.append(structure_set_roi)

    # ROI Contour Sequence
    roi_contour_sequence = Sequence()
    ds.ROIContourSequence = roi_contour_sequence

    # Loop over ROI contour sequences
    for ROI in range(1,numberOfROIs+1):

        # ROI Contour Sequence: ROI Contour 1
        roi_contour = Dataset()
        roi_contour.ROIDisplayColor = colors[ROI_color]

        # Contour Sequence
        contour_sequence = Sequence()
        roi_contour.ContourSequence = contour_sequence

        # Loop over slices in volume (ROI)
        for z in range(volume.shape[2]):

            # Should Contour Sequence be inside this loop?
            #contour_sequence = Sequence()
            #roi_contour.ContourSequence = contour_sequence

            # Loop over contour sequences in this slice
            numberOfContoursInThisSlice = len(allContourCoordinates[ROI-1][z])
            for c in range(numberOfContoursInThisSlice):

                currentCoordinates = allContourCoordinates[ROI-1][z][c]

                # Contour Sequence: Contour 1
                contour = Dataset()

                # Contour Image Sequence
                contour_image_sequence = Sequence()
                contour.ContourImageSequence = contour_image_sequence

                # Load the corresponding dicom file to get the SOPInstanceUID
                dstemp = pydicom.dcmread(input_dicom_path + "%s"%dicomFiles[z],stop_before_pixels=True)

                # Contour Image Sequence: Contour Image 1
                contour_image = Dataset()
                contour_image.ReferencedSOPClassUID = dstemp.SOPClassUID  # '1.2.840.10008.5.1.4.1.1.2'
                contour_image.ReferencedSOPInstanceUID = dstemp.SOPInstanceUID # '1.3.6.1.4.1.9590.100.1.2.76071554513024464020636223132290799275'
                contour_image_sequence.append(contour_image)

                contour.ContourGeometricType = 'CLOSED_PLANAR'
                contour.NumberOfContourPoints = len(currentCoordinates)
                contour.ContourData = currentCoordinates.tolist()
                contour_sequence.append(contour)

        roi_contour.ReferencedROINumber = ROI
        roi_contour_sequence.append(roi_contour)


    # RT ROI Observations Sequence
    rtroi_observations_sequence = Sequence()
    ds.RTROIObservationsSequence = rtroi_observations_sequence

    # Loop over ROI observations
    for ROI in range(1,numberOfROIs+1):
        # RT ROI Observations Sequence: RT ROI Observations 1
        rtroi_observations = Dataset()
        rtroi_observations.ObservationNumber = str(ROI)
        rtroi_observations.ReferencedROINumber = str(ROI)
        rtroi_observations.ROIObservationLabel = ''
        rtroi_observations.RTROIInterpretedType = ''
        rtroi_observations.ROIInterpreter = ''
        rtroi_observations_sequence.append(rtroi_observations)

    # Add RTSTRUCT specifics
    ds.Modality = 'RTSTRUCT' # So the GammaPlan software can recognize RTSTRUCT
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3' # So the GammaPlan software can recognize RTSTRUCT

    random_str_1 = "%0.8d" % random.randint(0,99999999)
    random_str_2 = "%0.8d" % random.randint(0,99999999)
    ds.SeriesInstanceUID = "1.2.826.0.1.3680043.2.1125."+random_str_1+".1"+random_str_2 # Just some random UID

    # Rename series for identification in GammaPlan
    ds.SeriesDescription = ROI_name

    # Save RTSTRUCT file
    ds.save_as(output_dicom_path + ROI_name + "_RTSTRUCT.dcm")

def get_parser():
    parser = argparse.ArgumentParser(description='Convert nifti segmentation volume to DICOM RTstruct')
    # Positional arguments.
    parser.add_argument("input_nifti", help="Path to input NIFTI image containing segmentation mask")
    parser.add_argument("input_dicom", help="Path to original DICOM files (used as template)")
    parser.add_argument("output_dicom", help="Path to output DICOM RTSTRUCT image")
    parser.add_argument("ROI_name", help="Name of this ROI")
    parser.add_argument("ROI_color", help="Color of this ROI")

    return parser.parse_args()

if __name__ == "__main__":
    p = get_parser()
    convert(p.input_nifti, p.input_dicom, p.output_dicom, p.ROI_name, p.ROI_color)
