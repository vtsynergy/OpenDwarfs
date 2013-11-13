#define CITIES 14

__kernel void search(
   __global  int *h,
   __global  int *city,
   __global int *result,
   __global int *traverse
){
   int i,j,counter,counter2,tmp,x,y,tmp2,tentative_score,tent_is_better,pos,start,end;
   int closedSet[CITIES];
   int openSet[CITIES];
   int f_score[CITIES];
   int g_score[CITIES];
   int visit[CITIES];
   int routine[CITIES];
   start = get_local_id(0);
   end = get_group_id(0);
   counter =1;
   counter2 =0;
   pos=-1;
   for(i=0;i<CITIES;i++){
      f_score[i]=0;
      g_score[i]=0;
      closedSet[i]=0;
      openSet[i]=0;
      visit[i]=0;
      routine[i]=0;
   }
  
   openSet[0] = start;
   g_score[0] = 0;
   f_score[0] = h[0];
   while(counter!=0){
      tmp = 0;
      tmp2 = f_score[openSet[0]];
      for(i=0;i<counter;i++){
          if(tmp2>f_score[openSet[i]]){
             tmp=i;
             tmp2=f_score[openSet[i]];
          }
      }

      x=openSet[tmp];
      if(x== end){
         result[get_local_id(0)+get_group_id(0)*CITIES]=g_score[x];
         tmp2=0;
         for(tmp=pos;x!=start;){
            routine[tmp2]=visit[x];
            x= visit[x];
            tmp2++;
         }

         tmp=tmp2-2;
         tmp2--;
         for(i=0;i<tmp2;i++){
            traverse[i+get_local_id(0)*CITIES+get_group_id(0)*CITIES*CITIES]=routine[tmp];
            tmp--;
         }
         traverse[i+get_local_id(0)*CITIES+get_group_id(0)*CITIES*CITIES]=-1;
         return;
      }

      closedSet[counter2]=x;
      counter2++;
      openSet[tmp]=openSet[counter-1];
      counter--;
      for(i=0;i<CITIES;i++){
          if(i==x){
             continue;
          }
          if(city[x*CITIES+i]==-1){
            continue;
          }
          y=i;
          tmp2=0;

          for(j=0;j<counter2;j++){
             if(y==closedSet[j]){
                tmp2=1;
             }
          }
          if(tmp2){
             continue;
          }
          tentative_score= g_score[x]+city[x*CITIES+y];
          tent_is_better =0;
          tmp2=0;

          for(j=0;j<counter;j++){
             if(y==openSet[j]){
                 tmp2=1;
             }
          }

          if(tmp2==0){
             openSet[counter]=i;
             counter++;
             tent_is_better=1;
          }else if(tentative_score<g_score[y]){
             tent_is_better=1;
          }
          else{
             tent_is_better=0;
          }
  
          if(tent_is_better ==1){
             visit[y]=x;
             pos++;
             g_score[y]=tentative_score;
             f_score[y]=g_score[y]+h[end*CITIES+y];
          }
      }

   }
   result[get_group_id(0)*CITIES+get_local_id(0)] =-1;
   return;
}