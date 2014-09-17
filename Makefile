pretty:
	find ./include \( -name "*.hpp" -o -name "*.h" -o -name "*.tcc" \) -type f -exec astyle -A2 '{}' \;
	find ./include -name "*.orig" -type f -exec rm '{}' \;
	find ./src \( -name "*.hpp" -o -name "*.h" -o -name "*.tcc" -o -name "*.cpp" \) -type f -exec astyle -A2 '{}' \;
	find ./src -name "*.orig" -type f -exec rm '{}' \;
