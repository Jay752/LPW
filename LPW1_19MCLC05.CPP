#include<stdio.h>
#include<conio.h>
#include<math.h>
void main()
{
	int i,J,k,u,NM,NJ,JE[12],KE[12],Dy[12],Rz[12],Im[12],Nl,Lcase,Rj[24],Lj,LJcase,nlj;
	float MOI[12],Elasticity[12],SM[5][5],SJ[12][12],l[12],Df[10],a[12],b[12],Dis[12],Am[24][24],ae[24],Aml[24],w,p,xl,m,temp=0,temp1=0,temp2=0,Jlv,Am1[24];

	FILE *fi,*fo;
	fi=fopen("Input.txt","r");
	fo=fopen("Output.txt","w");

	fscanf(fi,"%d %d",&NM,&NJ);
	fprintf(fo,"Number of member=%d\nNumber of Joint=%d",NM,NJ);

	fprintf(fo,"\n");
	for(i=1;i<=NJ;i++)
	{
		fscanf(fi,"%f",&Dis[i]);
		fprintf(fo,"\nFor Joint %d,Distance form 0th end is %0.2f",i,Dis[i]);
	}

	for(i=1;i<=NJ;i++)
	{
		fscanf(fi,"%d %d %d",&Dy[i],&Rz[i],&Rj[i]);
		fprintf(fo,"\nFor Joint %d,\tDisplacement[Y]=%d\tRotation[Z]=%d",i,Dy[i],Rz[i]);
	}

	for(i=1;i<=2*NJ;i++)
		for(J=1;J<=2*NJ;J++)
			SJ[i][J]=0.0;

	for(J=1;J<=2*NJ;J++)
		for(k=1;k<=2*NJ;k++)
			SM[J][k]=0.0;

	for(i=1;i<=NM;i++)
	{
		fscanf(fi,"%d %d %f %f",&JE[i],&KE[i],&MOI[i],&Elasticity[i]);
		fprintf(fo,"\n\n\t********For Member %d***********\n\nJEnd=%d\tKEnd=%d\tI=%0.2f\tE=%0.2f\t",i,JE[i],KE[i],MOI[i],Elasticity[i]);

		l[i]=(Dis[i+1]-Dis[i]);
		fprintf(fo,"\nLength of the member is %0.2f",l[i]);

		SM[1][1]=(12*Elasticity[i]*MOI[i]/pow(l[i],3));
		SM[1][2]=(6*Elasticity[i]*MOI[i]/pow(l[i],2));
		SM[1][3]=(-SM[1][1]);
		SM[1][4]=SM[1][2];
		SM[2][1]=SM[1][2];
		SM[2][2]=(4*Elasticity[i]*MOI[i]/l[i]);
		SM[2][3]=(-SM[2][1]);
		SM[2][4]=(SM[2][2]/2);
		SM[3][1]=(-SM[1][1]);
		SM[3][2]=(-SM[1][2]);
		SM[3][3]=(-SM[1][3]);
		SM[3][4]=(-SM[1][4]);
		SM[4][1]=SM[2][1];
		SM[4][2]=SM[2][4];
		SM[4][3]=SM[2][3];
		SM[4][4]=SM[2][2];

		fprintf(fo,"\n\n*****Stiffness matrix [K] [SM%d]=\n",i);
		for(J=1;J<=4;J++)
		{
			fprintf(fo,"\n");
			for(k=1;k<=4;k++)
			{
				fprintf(fo,"%0.2f\t",SM[J][k]);
			}
		}
		Im[1]=2*JE[i]-1;
		Im[2]=2*JE[i];
		Im[3]=2*KE[i]-1;
		Im[4]=2*KE[i];

		for(J=1;J<=4;J++)
			for(k=1;k<=4;k++)
				SJ[Im[J]][Im[k]]+=SM[J][k];
	}

	fprintf(fo,"\n\n******Joint stiffness matrix [Sj]=\n");
	for(i=1;i<=2*NJ;i++)
	{
		fprintf(fo,"\n");
		for(J=1;J<=2*NJ;J++)
		{
			fprintf(fo,"%0.2f\t",SJ[i][J]);
		}
	}

	for(i=1;i<2*NJ;i++)
	{
		ae[i]=0.0;
	}
	for(i=1;i<=NM;i++)
	{
		for(J=1;J<=4;J++)
		{
			Am[i][J]=0.0;
		}
		fscanf(fi,"%d",&Nl);
		for(J=1;J<=Nl;J++)
		{
			fscanf(fi,"%d",&Lcase);
			switch(Lcase)
			{
			case 1:
				fscanf(fi,"%f",&w);
				Am[i][1]+=w*l[i]/2;
				Am[i][2]+=w*l[i]*l[i]/12;
				Am[i][3]+=Am[i][1];
				Am[i][4]+=-Am[i][2];
				break;
			case 2:
				fscanf(fi,"%f %f",&p,&xl);
				Am[i][1]+=p*(l[i]-xl)/l[i];
				Am[i][2]+=p*(l[i]-xl)*(l[i]-xl)*xl/l[i]/l[i];
				Am[i][3]+=p*xl/l[i];
				Am[i][4]+=-p*(l[i]-xl)*xl*xl/l[i]/l[i];
				break;
			case 3:
				fscanf(fi,"%f %f",&m,&xl);
				Am[i][3]+=m*(l[i]-xl)/l[i]/l[i];
				Am[i][4]+=m*(l[i]-xl)*(l[i]-xl)*xl/l[i]/l[i]/l[i];
				Am[i][1]+=m*xl/l[i]/l[i];
				Am[i][2]+=-m*(l[i]-xl)*xl*xl/l[i]/l[i]/l[i];
				break;
			}
		}
		fprintf(fo,"\n\nAMl matrix for Member %d \n ",i);
		for(J=1;J<=4;J++)
		{
			fprintf(fo,"\nAMl[%d][%d]=%0.2f",i,J,Am[i][J]);
		}
	}

	for(i=1;i<=NM;i++)
	{
		Im[1]=2*JE[i]-1;
		Im[2]=2*JE[i];
		Im[3]=2*KE[i]-1;
		Im[4]=2*KE[i];

		for(J=1;J<=4;J++)
			ae[Im[J]]+=-Am[i][J];
	}

	fscanf(fi,"%d",&nlj);


	fprintf(fo,"\n");
	for(i=1;i<=nlj;i++)
	{
		fscanf(fi,"%d %d %f",&Lj,&LJcase,&Jlv);
		ae[2*Lj+LJcase-2]+=Jlv;
	}

	fprintf(fo,"\nCombined load vector [Ae]\n");
	for(i=1;i<=2*NJ;i++)
	{
		fprintf(fo,"\nAe[%d]=%0.2f",i,ae[i]);
	}

	for(i=1;i<=NJ;i++)
		{
		if(Dy[i]!=1)
		{
			SJ[2*Rj[i]-1][2*Rj[i]-1]+=1e8;
			ae[2*Rj[i]-1]+=Dy[i]*1e8;
		}
		if(Rz[i]!=1)
		{
			SJ[2*Rj[i]][2*Rj[i]]+=1e8;
			ae[2*Rj[i]]+=Rz[i]*1e8;
		}
	}

	for(u=1;u<=2*NJ;u++)
	{
		temp=SJ[u][u];
		for(J=1;J<=2*NJ;J++)
		{
			SJ[u][J]=SJ[u][J]/temp;
		}
		ae[u]=ae[u]/temp;
		for(i=1;i<=2*NJ;i++)
		{
			if(i==u) continue;
			temp1=SJ[i][u];
			for(J=1;J<=2*NJ;J++)
			{
				SJ[i][J]=SJ[i][J]-SJ[u][J]*temp1;
			}
			ae[i]=ae[i]-ae[u]*temp1;
		}
	}

	fprintf(fo,"\n\nDisplacement matrix [Delta]\n :");
	for(i=1;i<=2*NJ;i++)
	{
		fprintf(fo,"\ndJ[%d]=%0.2f",i,ae[i]);
	}


	for(i=1;i<=NM;i++)
	{
		Df[1]=ae[2*JE[i]-1];
		Df[2]=ae[2*JE[i]];
		Df[3]=ae[2*KE[i]-1];
		Df[4]=ae[2*KE[i]];

		for(J=1;J<=4;J++)
		{
			temp2=0;
			for(k=1;k<=4;k++)
			{
				temp2=temp2+SM[J][k]*Df[k];
				Am1[J]=temp2;
			}
		}
		fprintf(fo,"\n\n[AM] Member end action for member : %d \n",i);
		for(J=1;J<=4;J++)
			{
				Aml[J]=Am[i][J]+Am1[J];
				fprintf(fo,"\nAM[%d]=%0.2f",J,Aml[J]);
			}
		
	}			
	getch();
}