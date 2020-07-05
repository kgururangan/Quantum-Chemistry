function [] = print_intro(Nelec,nfzc,nfzv,Norb,nact_h,nact_p,run,input_file)


filename = strip(split(input_file,'.')); proj_name = filename{1};
	
	fprintf('\nProject Name: %s\n',proj_name)
	fprintf('---------------------------------\n')	
	fprintf('Working directory: %s\n\n',pwd')

	fprintf('System Information:\n')
	fprintf('---------------------\n')
	fprintf('   Number of electrons = %d\n',Nelec)
	fprintf('   Number of spatial orbitals = %d\n',Norb)
	fprintf('   Number of frozen core spatial orbitals = %d\n',nfzc)
	fprintf('   Number of frozen virtual spatial orbitals = %d\n',nfzv)
	fprintf('   Number of occupied spinorbitals = %d\n',Nelec-2*nfzc)
	fprintf('   Number of unoccupied spinorbitals = %d\n',2*Norb-Nelec-2*nfzv)
	fprintf('   Number of active occupied spatial orbitals = %d\n',nact_h)
	fprintf('   Number of active unoccupied spatial orbitals = %d\n',nact_p)

	if run.left_cc == 1
		fprintf('\nCalculation Type = %s + Left-CC\n',run.calc_type)
	else
		fprintf('\nCalculation Type = %s\n',run.calc_type)
	end

end
