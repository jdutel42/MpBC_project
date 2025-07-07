function remplace(strin,c1,c2)
{
var strout='';
var i=0;
while(i<strin.length) 
	{
	if (strin.charAt(i)==c1) strout=strout+c2;else strout=strout+strin.charAt(i);
	i++;
	}
return strout;
}

function testfloat(obj)
{
var valstr=obj.value;
valstr=remplace(valstr,',','');
valstr=remplace(valstr,' ','');
var val=0;
if (valstr!='') val=parseFloat(valstr);
if (isNaN(val)) val=0;
obj.value=val;
return val;
}

function testInt(obj)
{
var vallg=parseInt(obj.value,10);
if (isNaN(vallg)) vallg=0;
obj.value=vallg;
return vallg;
}


function majNOM(obj)
{
obj.value=obj.value.toUpperCase();
return obj.value;
}

function majPRENOM(obj)
{
var valpre=obj.value.toLowerCase();
var cprev=' ';
var ccur='';
var newpre='';
for (i=0;i<valpre.length;i++) {
	ccur=valpre.charAt(i);
	if ((cprev==' ')||(cprev=='-')||(cprev=='.')) ccur=ccur.toUpperCase();
	cprev=ccur;
	newpre=newpre+ccur;
	}
obj.value=newpre;
return newpre;
}


function ltrim(strin)
{
var strout=''+strin;
var i=0;
while((i<strout.length)&&((strout.charAt(i)==' ')||(strout.charAt(i)=='\r')||(strout.charAt(i)=='\n'))) i++;
if (i!=0) return strout.substr(i);
else return strout;
}
function rtrim(strin)
{
var strout=''+strin;
var i=strout.length-1;
while((i>=0)&&((strout.charAt(i)==' ')||(strout.charAt(i)=='\r')||(strout.charAt(i)=='\n'))) i--;
if (i!=strout.length-1) return strout.substr(0,i+1);
else return strout;
}

function btrim(strin)
{
return ltrim(rtrim(strin));
}

function getdate0(strin)
{

var strout=btrim(strin);
var lng=strout.length;
var i=0;

var jourstr='';
var cpt=0;while((i<lng) && isNaN(strout.charAt(i))) i++;
while ((i<lng) &&(cpt<2) && (!isNaN(strout.charAt(i)))) {jourstr=jourstr+strout.charAt(i);i++;cpt++;}
var jouri=parseInt(jourstr,10);if (isNaN(jouri)) jouri=0;
var moisstr='';
var cpt=0;while(isNaN(strout.charAt(i))) i++;
while ((i<lng) && (cpt<2) && (!isNaN(strout.charAt(i)))) {moisstr=moisstr+strout.charAt(i);i++;cpt++;}

var moisi=parseInt(moisstr,10);if (isNaN(moisi)) moisi=0;
var anneestr='';
var cpt=0;while(isNaN(strout.charAt(i))) i++;
while ((i<lng) && (cpt<4) && (!isNaN(strout.charAt(i)))) {anneestr=anneestr+strout.charAt(i);i++;cpt++;}
var anneei=parseInt(anneestr,10);if (isNaN(anneei)) anneei=0;
if (anneei<50) anneei=anneei+2000;
else if (anneei<100) anneei=anneei+1900;

if ((jouri<1)||(jouri>31)||(moisi<1)||(moisi>12)||(anneei<1900)||(anneei>2050)) return '';
else
	{
	var theDay=new Date(anneei,moisi-1,jouri);
//	theDay.setFullYear(anneei,moisi-1,jouri);
	return theDay;
	}
}

function getdate(strin)
{
var theDate=getdate0(strin);
if (theDate=='') return '';
else
	{
	var str=''+theDate.getDate()+'/'+(theDate.getMonth()+1)+'/'+theDate.getFullYear();
	return str;
	}
}

function openMax(theURL,winName) { //v2.0
  wid=window.screen.availWidth-10;
  hei=window.screen.availHeight-28;
window.open(theURL,winName,'top=0,left=0,status=yes,menubar=yes,width='+wid+'\',height='+hei+'\'');
}

function MM_openBrWindow(theURL,winName,features) { //v2.0
  window.open(theURL,winName,features);
}

function preloadImages() { 
  var d=document; if(d.images){ if(!d.p) d.p=new Array();
    var i,j=d.p.length,a=preloadImages.arguments; for(i=0; i<a.length; i++)
    if (a[i].indexOf("#")!=0){ d.p[j]=new Image; d.p[j++].src=a[i];}}
}

function findObj(n, d) {
  var p,i,x;
  if(document.getElementById)
	{
   	x=document.getElementById(n);
	return x;
	}
  else
	{
  	if(!d) d=document; if((p=n.indexOf("?"))>0&&parent.frames.length) {
    		d=parent.frames[n.substring(p+1)].document; n=n.substring(0,p);}
  	if(!(x=d[n])&&d.all) x=d.all[n]; for (i=0;!x&&i<d.forms.length;i++) x=d.forms[i][n];
  	for(i=0;!x&&d.layers&&i<d.layers.length;i++) x=findObj(n,d.layers[i].document); return x;
	}
}


function modifdiv(NSLayer,NSLayerCont,newval) { //v2.0
  var navappVer=navigator.appVersion;navappVer=navappVer[0];
  var NS = ((navigator.appName == 'Netscape') && (navappVer < 5));
  var obj;
  if (NS) {
	var objStr='document.layers[\''+NSLayerCont+'\'].document.layers[\''+NSLayer+'\']';
	var obj=eval(objStr);
	}
  else {obj=findObj(NSLayer);}
  if (NS && (obj!=null))
{
//obj.top=topVal;
obj.document.write('<p class=NSLAY>'+newval+'</p>');
obj.document.close();
obj.visibility='show';
}
  else if (obj != null) {
obj.style.visibility='visible';
obj.innerHTML=newval;
}
return true;
}


function MM_swapImgRestore() { //v2.0
  if (document.MM_swapImgData != null)
    for (var i=0; i<(document.MM_swapImgData.length-1); i+=2)
      document.MM_swapImgData[i].src = document.MM_swapImgData[i+1];
}

function MM_swapImage() { //v1.2
  var i,j=0,objStr,obj,swapArray=new Array,oldArray=document.MM_swapImgData;
  for (i=0; i < (MM_swapImage.arguments.length-2); i+=3) {
    objStr = MM_swapImage.arguments[(navigator.appName == 'Netscape')?i:i+1];
    if ((objStr.indexOf('document.layers[')==0 && document.layers==null) ||
        (objStr.indexOf('document.all[')   ==0 && document.all   ==null))
      objStr = 'document'+objStr.substring(objStr.lastIndexOf('.'),objStr.length);
    obj = eval(objStr);
    if (obj != null) {
      swapArray[j++] = obj;
      swapArray[j++] = (oldArray==null || oldArray[j-1]!=obj)?obj.src:oldArray[j];
      obj.src = MM_swapImage.arguments[i+2];
  } }
  document.MM_swapImgData = swapArray; //used for restore
}

function majCODE(obj,maxc)
{
var str=btrim(obj.value.toUpperCase());
var str0='';
var i=0;
while (i<str.length) {
	if (((str.charCodeAt(i)>=48)&&(str.charCodeAt(i)<=57))||(str.charCodeAt(i)==32)) break;
	else {
		str0=str0+str.charAt(i);
		i++;
		}
	}
if (i<str.length) {
	if ((i>0) && (str.charCodeAt(i)!=32)) str0=str0+' ';
	str0=str0+str.charAt(i);
	i++;
	while ((i<str.length) && (str.charAt(i)==' ')) i++;
	while (i<str.length) {
		//if (str.charAt(i)!=' ')	
		str0=str0+str.charAt(i);
		i++;
		}
	}
if (str0.length<=maxc) obj.value=str0;
return obj.value;
}


function clearForm() {
var	oForm=window.document.forms[0];
var frm_elements = oForm.elements;
for (i = 0; i < frm_elements.length; i++)
	{
    field_type = frm_elements[i].type.toLowerCase();
    switch (field_type)
    {
    case "text":
    case "password":
    case "textarea":
    case "hidden":
        frm_elements[i].value = "";
        break;
    case "radio":
    case "checkbox":
        if (frm_elements[i].checked)
        {
            frm_elements[i].checked = false;
        }
        break;
    case "select-one":
    case "select-multi":
        frm_elements[i].selectedIndex = -1;
        break;
    default:
        break;
    }
	}
}

