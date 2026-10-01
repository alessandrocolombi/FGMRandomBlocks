from pathlib import Path
import csv, math
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor
pdfmetrics.registerFont(TTFont('Arial','C:/Windows/Fonts/arial.ttf'))
pdfmetrics.registerFont(TTFont('ArialBold','C:/Windows/Fonts/arialbd.ttf'))
rows=list(csv.DictReader(open('tmp/pdfs/K_088_094_pmf.csv')))
s=next(csv.DictReader(open('tmp/pdfs/K_088_094_summary.csv')))
mu=float(s['mean']); sd=float(s['sd']); B=int(s['B'])
Path('output/pdf').mkdir(parents=True,exist_ok=True)
c=canvas.Canvas('output/pdf/PriorK_a_theta_0.88_b_theta_0.94.pdf',pagesize=(720,504))
c.setTitle('Prior marginale di K - a_theta=0.88, b_theta=0.94')
def text(x,y,t,size=10,font='Arial',color='#26354A'):
 c.setFillColor(HexColor(color)); c.setFont(font,size); c.drawString(x,y,t)
text(48,460,'Distribuzione prior di K',22,'ArialBold')
text(48,436,'σ ~ Beta(1, 1)     θ + σ | σ ~ Gamma(shape = 0.88, rate = 0.94)',11)
text(48,418,f'p = 40  •  {B:,} simulazioni Monte Carlo  •  9 sottogruppi fissati dagli eta'.replace(',',' '),10,color='#58697D')
for x,label,value in [(48,'MEDIA',f'{mu:.2f}'),(245,'DEVIAZIONE STANDARD',f'{sd:.2f}'),(480,'QUANTILI 2.5% / 97.5%',f"{int(float(s['q025']))} / {int(float(s['q975']))}")]:
 text(x,388,label,8,'ArialBold','#58697D');text(x,365,value,18,'ArialBold')
x0,y0,w,h=66,104,606,220
probs=[float(r['probability']) for r in rows if int(r['k'])>=9]
ymax=math.ceil(max(probs)*1.12/0.02)*0.02
step=w/32
for i in range(round(ymax/0.02)+1):
 val=i*.02; yy=y0+h*val/ymax
 c.setStrokeColor(HexColor('#DFE5EC')); c.setLineWidth(.5); c.line(x0,yy,x0+w,yy)
 text(30,yy-3,f'{100*val:.0f}%',9,color='#58697D')
for j,p in enumerate(probs):
 c.setFillColor(HexColor('#327FAD'));c.rect(x0+j*step+1.5,y0,step-3,h*p/ymax,fill=1,stroke=0)
for k in [9,10,15,20,25,30,35,40]:
 xx=x0+(k-9+.5)*step
 text(xx-5,y0-17,str(k),9)
c.setStrokeColor(HexColor('#D17B25'));c.setLineWidth(1.5);c.setDash(4,3)
xx=x0+(mu-9+.5)*step;c.line(xx,y0,xx,y0+h)
c.setDash();text(478,337,'Linea tratteggiata: media di K',9,color='#A15A19')
text(x0,337,'Probabilità stimata P(K = k)',10,'ArialBold')
text(269,65,'Numero totale di cluster K',11)
text(48,39,'Eta = 1 nelle posizioni 4, 6, 9, 13, 18, 22, 28, 33, 40; zero altrove.',9,color='#58697D')
text(48,24,'Barre = frequenze relative; nessuna lisciatura. Supporto: K = 9, …, 40. Seme: 20261001.',9,color='#58697D')
c.showPage();c.save()
