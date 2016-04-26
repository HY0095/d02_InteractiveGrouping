  case when A1 <= 20.0 then 0.329162622592
       when A1 > 20.0 and A1 <= 30.0 then 0.121240574729
       when A1 > 30.0 and A1 <= 59.0 then -0.0978037420699
       when A1 > 59.0 then -0.326063381827
       else 0.0000 end as A1_woe
 ,case when A127 <= 3.0 then 0.383811346501
       when A127 > 3.0 and A127 <= 7.0 then 0.0203666791749
       when A127 > 7.0 and A127 <= 10.0 then -0.382455418218
       when A127 > 10.0 then -0.768147068415
       else 0.0000 end as A127_woe
 ,case when A86 <= 2000.0 then -0.000750250963242
       when A86 > 2000.0 and A86 <= 17000.0 then 0.0783652222197
       when A86 > 17000.0 and A86 <= 24000.0 then -0.174954196963
       when A86 > 24000.0 then -0.40035946366
       else 0.0000 end as A86_woe
 ,case when A146 <= 1.0 then 0.0644589095937
       when A146 > 1.0 and A146 <= 2.0 then -0.150868533057
       when A146 > 2.0 then -0.568927014833
       else 0.0000 end as A146_woe
