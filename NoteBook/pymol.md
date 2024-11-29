# temp

（1）pymol中选中残基编号为负数的单个残基。需要在编号前面加\。   
```python
PyMOL>select chain I and (resi \-61)
 Selector: selection "sele" defined with 19 atoms.
PyMOL>select chain I and (resi \1)
 Selector: selection "sele" defined with 22 atoms.
```