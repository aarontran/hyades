media/column.mp4: flds.all.hdf5
	./column.py
	./column_makemov.sh

flds.all.hdf5: ../output collate.py
	./collate.py ../output flds

clean:
	rm -f flds.all.hdf5
	rm -f media/column.mp4
