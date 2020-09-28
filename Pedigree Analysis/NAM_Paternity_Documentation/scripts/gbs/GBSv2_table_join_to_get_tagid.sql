--From Liang Gao

SELECT  	sp.snpid,
			sp.chromosome, 
			sp.position,
			a.snpid,
			a.alleleid,
			ta.tagid
FROM snpposition sp 
JOIN allele a
		on sp.snpid = a.snpid
JOIN tagallele ta
		on a.alleleid=ta.alleleid
order by sp.chromosome, sp.position;
