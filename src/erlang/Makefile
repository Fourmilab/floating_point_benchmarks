
fbench.beam:	fbench.erl
	erlc fbench.erl
	
clean:
	rm -f fbench.beam erl_crash.dump
	
time:	fbench.beam
	time -p erl -noinput fbench.beam -s fbench main -s init stop >results.txt
	diff results.txt correct_results.txt ; true

#	The following targets will only work if you have
#	an Erlang/OTP system installed which supports HiPE
#	(High Performance Erlang) native code compilation.
#	HIPE should be set to the binaries directory of the
#	installation.

HiPE = ${HOME}/linuxtools/otp_src_17.0/bin

native:
	$(HiPE)/erlc +native '+{hipe, [o3]}' fbench.erl

native_time:
	time -p $(HiPE)/erl -noinput fbench.beam -s fbench main -s init stop >results.txt
	diff results.txt correct_results.txt ; true

hipe:
	$(HiPE)/erlc fbench.erl
	
hipe_time:
	time -p $(HiPE)/erl -noinput fbench.beam -s fbench main -s init stop >results.txt
	diff results.txt correct_results.txt ; true

		
