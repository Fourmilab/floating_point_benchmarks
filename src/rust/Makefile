
duh:
	@echo "What'll it be--clean, debug, or time?"

clean:
	( rm -f rm RustFB_results.txt ; cd fbench ; cargo clean )

debug:
	( cd fbench ; CARGO_INCREMENTAL=0 cargo build )
	echo "Press return to begin benchmark: Stop the timer: " >RustFB_results.txt
	fbench/target/debug/fbench | head --lines=-1 >>RustFB_results.txt
	diff RustFB_results.txt ../pascal/results.txt
	rm RustFB_results.txt

time:
	( cd fbench ; CARGO_INCREMENTAL=0 cargo build --release )
	( time -p fbench/target/release/fbench )
