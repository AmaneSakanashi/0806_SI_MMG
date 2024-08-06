# SI-CMAES-MMG-esso
System Identification of MMG model for Esso Osaka using CMA-ES

Archive purpose only.
If you want to use CMA-ES for SI, use [Abkowtiz model version](https://github.com/NAOE-5thLab/SI-CMA-abkowitz-esso)


##　ビルドする時
    cmake .. -G "Unix Makefiles" -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Debug
    make  
## リビルドする時
    cd ../
    rm -rf ./build/*
    cd ./build
    cmake .. -G "Unix Makefiles" -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release
    make
##　実行
./main.out ../inputfiles/InitCmaes57.txt  ../inputfiles/InitProblem57.txt ../inputfiles/random_train_file_list.txt 
