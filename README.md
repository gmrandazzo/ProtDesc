# ProtDesc - Protein Descriptors

Read a fasta file and generate a fingerprint based on triad of aminoacids

![ScreenShot](https://raw.githubusercontent.com/gmrandazzo/ProtDesc/main/triad_example.png)


## Usage

- makeprotdb.py: Create a database of triads
- protodes.py: Build descriptors based on triads features

## Fast database creation

create a script mkdb.x 

mkdb.x
```
python3 ~/Nextcloud/Software/ProtDesc/src/makeprotdb.py ${1} ${1}.db
```

find .  -name "*.fasta" | xargs -n 1 -P 12 -I {} bash mkdb.x {}

## Fast descriptor calculation

create a script mkdesc.x 


```
python3 /home/marco/Nextcloud/Software/ProtDesc/src/protdesc.py ${1} /home/marco/Nextcloud/Software/ProtDesc/src/db ${1}.desc.csv
```

find .  -name "*.fasta" | xargs -n 1 -P 12 -I {} bash mkdesc.x {}


## Info

**Author: Giuseppe Marco Randazzo <br/>
Mantainer: Giuseppe Marco Randazzo, gmrandazzo@gmail.com <br/>**

License
-------

ProtDesc is distributed under GPLv3 license.
For more details please read the file "LICENSE" or go to "http://www.gnu.org/licenses/gpl-3.0.html"



