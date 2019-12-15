

TGraphErrors* Residuals(TGraphErrors* data, TF1* f, bool fit = false, bool draw = false){
    TGraphErrors* res = new TGraphErrors(*data);
    if(fit)
        data->Fit(f);
    
    for(unsigned int i = 0; i < data->GetN(); i++)
        *(res->GetY() + i) -= f->Eval(*(res->GetX() + i));
    
    if(draw){
        res->Draw();
        TF1* zero = new TF1("zero","[0]*x",res->GetXaxis()->GetXmin(),res->GetXaxis()->GetXmax());
        zero->SetParameter(0,0);
        zero->Draw("SAME");
    }
    
    return res;
}
