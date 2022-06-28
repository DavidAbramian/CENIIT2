
# subj=$1
str1="z5_c0"

# series_num=1

# for subj in "17904" "18582" "18975" "19015" "19849" ; do
for subj in "19849" ; do

  series_num=1

  for contourLevel in "0.5" "0.4" "0.6" ; do

    for str2 in "6mm" "6mm_m" "6mm_m2" "6mm_m4" "4mm" "4mm_m" "4mm_m2" "4mm_m4" ; do

      series_num=$(( series_num + 1 ))

      echo $subj $str1 $str2 $series_num $contourLevel
      python convert_to_rtstruct_multi_ROI_multi_file_2.py $subj $str1 $str2 $series_num -l $contourLevel
      echo

    done

  done

done