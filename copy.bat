@echo off

set simulator_folder=..\ogs_kb1
set contraflow_folder=%simulator_folder%\Libs\Contraflow

md %contraflow_folder%

copy build\src\Release\contraflow.lib %contraflow_folder%
copy src\*.h %contraflow_folder%
copy src\matrix_stru3\*.h %contraflow_folder%


